//! Winnow lexical combinators plus a bounded Pratt parser. No error recovery on run.
use super::ast::*;
use anyhow::{bail, Context, Result};
use winnow::{combinator::alt, token::take_while, Parser};

const MAX_BYTES: usize = 256 * 1024;
const MAX_TOKENS: usize = 16_384;
const MAX_DEPTH: usize = 64;

#[derive(Debug)]
struct Token {
    text: String,
    at: usize,
}

fn lexeme<'a>(input: &mut &'a str) -> winnow::Result<&'a str> {
    alt((
        alt(("|>", "..", ">>", "~>", "==", "!=", "<=", ">=", "->")),
        take_while(1.., |c: char| c.is_ascii_alphanumeric() || c == '_'),
        take_while(1, |c: char| "@{}[]():,;.=!&|+-*/<>".contains(c)),
    ))
    .parse_next(input)
}

fn lex(source: &str) -> Result<Vec<Token>> {
    if source.len() > MAX_BYTES {
        bail!("GQ source exceeds {MAX_BYTES} bytes");
    }
    let mut rest = source;
    let mut tokens = Vec::new();
    while !rest.is_empty() {
        rest = rest.trim_start();
        if rest.is_empty() {
            break;
        }
        if rest.starts_with('#') {
            rest = rest.find('\n').map_or("", |n| &rest[n..]);
            continue;
        }
        let at = source.len() - rest.len();
        let text = if rest.starts_with('"') {
            let mut escaped = false;
            let mut end = None;
            for (i, ch) in rest.char_indices().skip(1) {
                if ch == '"' && !escaped {
                    end = Some(i + 1);
                    break;
                }
                escaped = ch == '\\' && !escaped;
            }
            let end = end.with_context(|| format!("byte {at}: unterminated string"))?;
            let text = rest[..end].to_owned();
            let _: String = serde_json::from_str(&text).context("invalid string escape")?;
            rest = &rest[end..];
            text
        } else {
            lexeme
                .parse_next(&mut rest)
                .map_err(|_| anyhow::anyhow!("byte {at}: unexpected character"))?
                .to_owned()
        };
        tokens.push(Token { text, at });
        if tokens.len() > MAX_TOKENS {
            bail!("GQ token budget exceeds {MAX_TOKENS}");
        }
    }
    tokens.push(Token {
        text: String::new(),
        at: source.len(),
    });
    Ok(tokens)
}

struct Cursor {
    tokens: Vec<Token>,
    i: usize,
    depth: usize,
}
impl Cursor {
    fn peek(&self) -> &str {
        &self.tokens[self.i].text
    }
    fn at(&self) -> usize {
        self.tokens[self.i].at
    }
    fn bump(&mut self) -> String {
        let s = self.peek().to_owned();
        if self.i + 1 < self.tokens.len() {
            self.i += 1;
        }
        s
    }
    fn eat(&mut self, s: &str) -> bool {
        if self.peek() == s {
            self.bump();
            true
        } else {
            false
        }
    }
    fn need(&mut self, s: &str) -> Result<()> {
        if !self.eat(s) {
            bail!(
                "byte {}: expected {s:?}, found {:?}",
                self.at(),
                self.peek()
            );
        }
        Ok(())
    }
    fn name(&mut self) -> Result<String> {
        if !self
            .peek()
            .starts_with(|c: char| c.is_ascii_alphabetic() || c == '_')
            || !self
                .peek()
                .chars()
                .all(|c| c.is_ascii_alphanumeric() || c == '_')
        {
            bail!("byte {}: expected identifier", self.at());
        }
        Ok(self.bump())
    }
    fn unsigned(&mut self) -> Result<u32> {
        self.bump()
            .parse()
            .with_context(|| format!("byte {}: expected u32 coordinate", self.at()))
    }
    fn type_name(&mut self) -> Result<String> {
        let mut name = self.name()?;
        if self.eat("<") {
            name.push('<');
            let mut count = 0;
            loop {
                count += 1;
                if count > 4 {
                    bail!("type has too many qualifiers");
                }
                name.push_str(&self.name()?);
                if self.eat(">") {
                    name.push('>');
                    break;
                }
                self.need(",")?;
                name.push(',');
            }
        }
        Ok(name)
    }
    fn fields(&mut self, assigned: bool) -> Result<Vec<(String, Expr)>> {
        self.need("{")?;
        let mut fields = Vec::new();
        while !self.eat("}") {
            let expr = self.expr(0)?;
            let (name, value) = if self.eat("=") {
                let Kind::Name(name) = expr.kind else {
                    bail!("byte {}: field name must be identifier", expr.at);
                };
                (name, self.expr(0)?)
            } else {
                if assigned {
                    bail!("byte {}: expected named field = expression", expr.at);
                }
                (field_name(&expr)?, expr)
            };
            if fields.iter().any(|(n, _)| n == &name) {
                bail!("duplicate field {name}");
            }
            fields.push((name, value));
            if self.peek() != "}" {
                self.need(",")?;
            }
        }
        if fields.is_empty() {
            bail!("empty field list");
        }
        Ok(fields)
    }
    fn args(&mut self) -> Result<Vec<(Option<String>, Expr)>> {
        self.need("(")?;
        let mut args = Vec::new();
        while !self.eat(")") {
            let name = if self.tokens.get(self.i + 1).is_some_and(|t| t.text == ":") {
                let name = self.name()?;
                self.need(":")?;
                Some(name)
            } else {
                None
            };
            args.push((name, self.expr(0)?));
            if self.peek() != ")" {
                self.need(",")?;
            }
        }
        Ok(args)
    }
    fn expr(&mut self, minimum: u8) -> Result<Expr> {
        self.depth += 1;
        if self.depth > MAX_DEPTH {
            bail!("byte {}: expression nesting exceeds {MAX_DEPTH}", self.at());
        }
        let at = self.at();
        let token = self.bump();
        let kind = match token.as_str() {
            "+" | "-" if self.peek() == ")" => Kind::String(token),
            "!" | "-" => Kind::Unary(token, Box::new(self.expr(11)?)),
            "(" => {
                let e = self.expr(0)?;
                self.need(")")?;
                e.kind
            }
            "@" => Kind::Resource(self.name()?),
            "true" => Kind::Bool(true),
            "false" => Kind::Bool(false),
            "null" => Kind::Null,
            "unknown" => Kind::Unknown,
            "any" | "all" if self.peek() != "|>" && !self.peek().is_empty() => {
                let domain = self.name()?;
                if domain == "placement" {
                    bail!("byte {at}: placement is reserved; choose unique, stored (diagnostic), or multimap {{ ... alternative {{ ... }} }}");
                }
                if !["unique", "stored", "multimap", "alternative", "record"]
                    .contains(&domain.as_str())
                {
                    bail!("byte {at}: unknown quantifier domain {domain}");
                }
                self.need("{")?;
                let body = Box::new(self.expr(0)?);
                self.need("}")?;
                Kind::Quant {
                    all: token == "all",
                    domain,
                    body,
                }
            }
            "path" if self.peek() == "(" => {
                self.need("(")?;
                let order = if self.eat("order") {
                    self.need(":")?;
                    self.name()?
                } else {
                    self.bump()
                };
                if !["genomic", "+", "-"].contains(&order.as_str()) {
                    bail!("byte {at}: path order must be genomic, +, or -");
                }
                self.need(")")?;
                self.need("{")?;
                let body = Box::new(self.expr(0)?);
                self.need("}")?;
                Kind::Path { order, body }
            }
            "g" | "j" if self.eat("[") => {
                let mut chrom = self.bump();
                if chrom.starts_with('"') {
                    chrom = serde_json::from_str(&chrom)?;
                }
                self.need(":")?;
                let start = self.unsigned()?;
                self.need("..")?;
                let end = self.unsigned()?;
                if start >= end {
                    bail!("byte {at}: coordinates must satisfy start < end (0-based, half-open)");
                }
                let reverse = if self.eat(":") {
                    match self.bump().as_str() {
                        "+" => Some(false),
                        "-" => Some(true),
                        _ => bail!("byte {at}: strand must be + or -"),
                    }
                } else {
                    None
                };
                self.need("]")?;
                Kind::Locus {
                    junction: token == "j",
                    chrom,
                    start,
                    end,
                    reverse,
                }
            }
            "[" => {
                let mut list = Vec::new();
                while !self.eat("]") {
                    list.push(self.expr(0)?);
                    if self.peek() != "]" {
                        self.need(",")?;
                    }
                }
                Kind::List(list)
            }
            _ if token.starts_with('"') => Kind::String(serde_json::from_str(&token)?),
            _ if token.starts_with(|c: char| c.is_ascii_digit()) => {
                if let Some(s) = token.strip_suffix("u64") {
                    Kind::Count(s.parse().context("invalid UInt64 count literal")?)
                } else if let Some(s) = token.strip_suffix("bp") {
                    Kind::Bases(s.parse().context("invalid base quantity")?)
                } else {
                    let mut number = token;
                    if self.eat(".") {
                        number.push('.');
                        number.push_str(&self.bump());
                    }
                    let value: f64 = number.parse().context("invalid numeric literal")?;
                    if !value.is_finite() {
                        bail!("byte {at}: non-finite number");
                    }
                    Kind::Number(value)
                }
            }
            _ if token.starts_with(|c: char| c.is_ascii_alphabetic() || c == '_') => {
                if self.peek() == "(" {
                    Kind::Call(token, self.args()?)
                } else {
                    Kind::Name(token)
                }
            }
            _ => bail!("byte {at}: expected expression, found {token:?}"),
        };
        let mut lhs = Expr { at, kind };
        let mut operators = 0;
        loop {
            operators += 1;
            if operators > 128 {
                bail!("byte {at}: flat expression exceeds 128 operators");
            }
            if self.eat(".") {
                lhs = Expr {
                    at,
                    kind: Kind::Field(Box::new(lhs), self.name()?),
                };
                continue;
            }
            let op = self.peek();
            let precedence = match op {
                "|" => 1,
                "&" => 3,
                "==" | "!=" | "<" | ">" | "<=" | ">=" | "in" | "not" => 5,
                "+" | "-" => 7,
                "*" | "/" => 9,
                ">>" | "~>" => 13,
                _ => break,
            };
            if precedence < minimum {
                break;
            }
            let mut op = self.bump();
            if op == "not" {
                self.need("in")?;
                op = "not in".into();
            }
            let rhs = self.expr(precedence + 1)?;
            lhs = Expr {
                at,
                kind: Kind::Binary(op, Box::new(lhs), Box::new(rhs)),
            };
        }
        self.depth -= 1;
        Ok(lhs)
    }
}

pub fn field_name(e: &Expr) -> Result<String> {
    match &e.kind {
        Kind::Name(n) => Ok(n.clone()),
        Kind::Field(base, field) => Ok(format!("{}.{}", field_name(base)?, field)),
        _ => bail!("byte {}: computed field requires name = expression", e.at),
    }
}

pub fn parse(source: &str) -> Result<Document> {
    let mut p = Cursor {
        tokens: lex(source)?,
        i: 0,
        depth: 0,
    };
    p.need("header")?;
    let header = p.fields(true)?;
    let mut bindings = Vec::new();
    let mut functions = Vec::new();
    while p.peek() != "from" {
        if p.eat(";") {
            continue;
        }
        if p.eat("let") {
            let name = p.name()?;
            p.need("=")?;
            bindings.push((name, p.expr(0)?));
        } else {
            let exported = p.eat("export");
            p.need("fn")?;
            let name = p.name()?;
            p.need("(")?;
            let mut parameters = Vec::new();
            while !p.eat(")") {
                let n = p.name()?;
                let ty = if p.eat(":") {
                    Some(p.type_name()?)
                } else {
                    None
                };
                parameters.push((n, ty));
                if p.peek() != ")" {
                    p.need(",")?;
                }
            }
            let returns = if p.eat("->") {
                Some(p.type_name()?)
            } else {
                None
            };
            if exported && (returns.is_none() || parameters.iter().any(|(_, t)| t.is_none())) {
                bail!("export fn requires parameter and return types");
            }
            p.need("=")?;
            functions.push(Function {
                name,
                parameters,
                returns,
                exported,
                body: p.expr(0)?,
            });
        }
    }
    p.need("from")?;
    let source = p.expr(0)?;
    let mut stages = Vec::new();
    while p.eat("|>") {
        stages.push(match p.bump().as_str() {
            "support" => {
                p.need("(")?;
                p.need("unit")?;
                p.need(":")?;
                let unit = p.name()?;
                let by = if p.eat(",") {
                    p.need("by")?;
                    p.need(":")?;
                    p.fields(false)?
                } else {
                    Vec::new()
                };
                p.need(")")?;
                Stage::Support(unit, by)
            }
            "within" => Stage::Within(p.expr(0)?),
            "where" => Stage::Where(p.expr(0)?),
            "derive" => Stage::Derive(p.fields(true)?),
            "tally" => {
                let fields = p.fields(false)?;
                let by = if p.eat("by") {
                    p.fields(false)?
                } else {
                    Vec::new()
                };
                Stage::Tally(fields, by)
            }
            "summarize" => {
                let fields = p.fields(true)?;
                let by = if p.eat("by") {
                    p.fields(false)?
                } else {
                    Vec::new()
                };
                Stage::Summarize(fields, by)
            }
            "select" => Stage::Select(p.fields(false)?),
            "sort" => Stage::Sort(p.fields(false)?),
            "take" => Stage::Take(p.unsigned()? as usize),
            name => bail!("byte {}: unknown pipeline stage {name:?}", p.at()),
        });
    }
    p.eat(";");
    p.need("")?;
    let document = Document {
        header,
        bindings,
        functions,
        source,
        stages,
    };
    check_depth(&document)?;
    Ok(document)
}

fn check_depth(d: &Document) -> Result<()> {
    let mut stack = Vec::new();
    stack.push((&d.source, 0));
    stack.extend(d.header.iter().chain(&d.bindings).map(|(_, e)| (e, 0)));
    stack.extend(d.functions.iter().map(|f| (&f.body, 0)));
    for stage in &d.stages {
        match stage {
            Stage::Within(e) | Stage::Where(e) => stack.push((e, 0)),
            Stage::Derive(fields)
            | Stage::Select(fields)
            | Stage::Sort(fields)
            | Stage::Support(_, fields) => stack.extend(fields.iter().map(|(_, e)| (e, 0))),
            Stage::Tally(fields, by) | Stage::Summarize(fields, by) => {
                stack.extend(fields.iter().chain(by).map(|(_, e)| (e, 0)))
            }
            Stage::Take(_) => {}
        }
    }
    while let Some((e, depth)) = stack.pop() {
        if depth > 128 {
            bail!("byte {}: AST depth exceeds 128", e.at);
        }
        match &e.kind {
            Kind::Unary(_, a)
            | Kind::Field(a, _)
            | Kind::Quant { body: a, .. }
            | Kind::Path { body: a, .. } => stack.push((a, depth + 1)),
            Kind::Binary(_, a, b) => {
                stack.push((a, depth + 1));
                stack.push((b, depth + 1));
            }
            Kind::List(items) => stack.extend(items.iter().map(|e| (e, depth + 1))),
            Kind::Call(_, args) => stack.extend(args.iter().map(|(_, e)| (e, depth + 1))),
            _ => {}
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn deterministic_mutations_and_large_inputs_fail_without_panics() {
        let seeds = [include_str!("../../../../examples/gq/cd45.gq"),
            "header {gq=1,assembly=\"test\"} from @x.records |> within all |> summarize {n=count()}"];
        let alphabet: Vec<char> = "@{}[]():,;.=!&|+-*/<>\"\\#\n09é🧬".chars().collect();
        let mut state = 0x97531_u64;
        for i in 0..4096 {
            let mut next = || {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                state as usize
            };
            let mut chars: Vec<char> = seeds[i % seeds.len()].chars().collect();
            for _ in 0..1 + next() % 12 {
                let index = next() % chars.len().max(1);
                match next() % 3 {
                    0 if !chars.is_empty() => {
                        chars.remove(index);
                    }
                    1 if !chars.is_empty() => chars[index] = alphabet[next() % alphabet.len()],
                    _ => chars.insert(index, alphabet[next() % alphabet.len()]),
                }
            }
            let source: String = chars.into_iter().collect();
            if let Ok(doc) = parse(&source) {
                if let Ok(compiler) =
                    crate::gq::check::Compiler::new(&doc, std::collections::BTreeMap::new())
                {
                    let _ = compiler.compile();
                }
            }
        }
        for expression in ["true | ".repeat(7000) + "true", "!".repeat(1000) + "true"] {
            assert!(parse(&format!(
                "header {{gq=1}} from @x.records |> where {expression}"
            ))
            .is_err());
        }
        assert!(parse(&" ".repeat(MAX_BYTES + 1)).is_err());
    }
    #[test]
    #[ignore = "manual parser-only benchmark; run in release mode"]
    fn parser_benchmark() {
        let source = include_str!("../../../../examples/gq/cd45.gq");
        let iterations = 10_000;
        let start = std::time::Instant::now();
        for _ in 0..iterations {
            std::hint::black_box(parse(std::hint::black_box(source)).unwrap());
        }
        println!(
            "GQ parser: {} ns/document; {iterations} iterations; {} bytes/document",
            start.elapsed().as_nanos() / iterations,
            source.len()
        );
    }
    #[test]
    fn examples_and_precedence() {
        let q = parse("header {gq=1} let A=j[1:10..20:-] from @x.records |> within g[1:1..90] |> derive {a=any unique{A}, b=all stored{!A|A&unknown}} |> tally {a,b} by {sample,cell.type}").unwrap();
        assert_eq!(q.stages.len(), 3);
        assert_eq!(q.bindings.len(), 1);
    }
    #[test]
    fn invalid_documents_fail_closed() {
        for s in [
            "header {gq=1} from @x.records |> where any placement{true}",
            "header {gq=1} from @x.records garbage",
            "header {gq=1} from @x.records |> within g[1:20..10]",
            "header {gq=1} from @x.records |> where \"é",
        ] {
            assert!(parse(s).is_err(), "{s}");
        }
        let deep = format!(
            "header {{gq=1}} from @x.records |> where {}true{}",
            "(".repeat(100),
            ")".repeat(100)
        );
        assert!(parse(&deep).is_err());
    }
    #[test]
    fn comments_and_unicode_strings() {
        assert!(parse("header {gq=1} # hi\nfrom @x.cells |> where cell.type == \"neuröne\" |> summarize {n=count()}").is_ok());
    }
}
