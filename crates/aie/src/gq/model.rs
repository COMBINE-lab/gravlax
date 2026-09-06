//! Typed logical expressions shared by planning and native execution.
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Truth {
    False,
    Unknown,
    True,
}
impl Truth {
    pub fn not(self) -> Self {
        match self {
            Self::True => Self::False,
            Self::False => Self::True,
            Self::Unknown => self,
        }
    }
    pub fn and(self, b: Self) -> Self {
        self.min(b)
    }
    pub fn or(self, b: Self) -> Self {
        self.max(b)
    }
    pub fn from_bool(b: bool) -> Self {
        if b {
            Self::True
        } else {
            Self::False
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            Self::False => "false",
            Self::True => "true",
            Self::Unknown => "unknown",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Type {
    Truth,
    String,
    Number,
    Count,
    Bases,
    Null,
    Region,
    Pattern,
    PatternSet,
    GeneFeature,
    TranscriptFeature,
    Bounds,
    List(Box<Type>),
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Region {
    pub chrom: String,
    pub intervals: Vec<(u32, u32)>,
    pub reverse: Option<bool>,
}
impl Region {
    pub fn contains(&self, p: u32) -> bool {
        self.intervals.iter().any(|&(a, b)| a <= p && p < b)
    }
    pub fn union(&mut self) {
        self.intervals.sort_unstable();
        let mut out: Vec<(u32, u32)> = Vec::new();
        for &(a, b) in &self.intervals {
            if let Some(last) = out.last_mut().filter(|last| last.1 >= a) {
                last.1 = last.1.max(b);
            } else {
                out.push((a, b));
            }
        }
        self.intervals = out;
    }
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Pattern {
    pub chrom: String,
    pub reverse: Option<bool>,
    pub junctions: Vec<(u32, u32)>,
    pub subsequence: bool,
    pub left: u32,
    pub right: u32,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub enum Value {
    Truth(Truth),
    String(String),
    Number(f64),
    Count(u64),
    Bases(u32),
    Null,
    Region(Region),
    Pattern(Pattern),
    Bounds { lower: u64, upper: u64 },
    List(Vec<Value>),
    PatternSet(Vec<Pattern>),
    Feature(Box<Feature>),
}
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Feature {
    pub transcript: bool,
    pub span: Region,
    pub exons: Region,
    pub junctions: Vec<Pattern>,
    pub junction_path: Option<Pattern>,
}
impl Value {
    pub fn truth(&self) -> Result<Truth> {
        match self {
            Self::Truth(v) => Ok(*v),
            Self::Null => Ok(Truth::Unknown),
            _ => bail!("expected Truth"),
        }
    }
    pub fn json(&self) -> serde_json::Value {
        use serde_json::json;
        match self {
            Self::Truth(v) => json!(v.name()),
            Self::String(v) => json!(v),
            Self::Number(v) => json!(v),
            Self::Count(v) => json!(v),
            Self::Bases(v) => json!(v),
            Self::Null => serde_json::Value::Null,
            Self::Bounds { lower, upper } => {
                json!({"lower":lower,"upper":upper,"status":if lower==upper{"exact"}else{"bounded"}})
            }
            Self::List(v) => v.iter().map(Self::json).collect(),
            _ => json!(self),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Unit {
    Record,
    Class,
    Cell,
    Alignment,
    Signature,
}
#[derive(Clone, Debug, Serialize)]
pub struct Node {
    pub at: usize,
    pub ty: Type,
    pub op: Op,
}
#[derive(Clone, Debug, Serialize)]
pub enum Op {
    Constant(Value),
    Field(String),
    Unary(String, Box<Node>),
    Binary(String, Box<Node>, Box<Node>),
    Quant {
        all: bool,
        domain: String,
        body: Box<Node>,
    },
    Geometry {
        kind: String,
        region: Region,
        minimum: u32,
        total: bool,
    },
    Match(Pattern),
    TxStrand(bool),
    Terminal(Region),
    Nonempty(String),
    Call(String, Vec<Node>),
    Support(Box<Node>),
}
#[derive(Clone, Debug, Serialize)]
pub enum Stage {
    Where(Node),
    Derive(Vec<(String, Node)>),
    Tally(Vec<(String, Node)>, Vec<(String, Node)>),
    Summarize(Vec<(String, Node)>, Vec<(String, Node)>),
    Select(Vec<(String, Node)>),
    Sort(Vec<(String, Node)>),
    Take(usize),
}
#[derive(Clone, Debug, Serialize)]
pub struct Plan {
    pub enumeration: bool,
    pub source: String,
    pub unit: Unit,
    pub assembly: String,
    pub library: Option<String>,
    pub within: Option<Value>,
    pub stages: Vec<Stage>,
    pub terminal: bool,
    pub diagnostic_stored: bool,
    pub fields: BTreeMap<String, Type>,
}

pub fn scalar_binary(op: &str, a: Value, b: Value) -> Result<Value> {
    if op == "&" {
        return Ok(Value::Truth(a.truth()?.and(b.truth()?)));
    }
    if op == "|" {
        return Ok(Value::Truth(a.truth()?.or(b.truth()?)));
    }
    if a == Value::Null || b == Value::Null {
        return Ok(if ["+", "-", "*", "/"].contains(&op) {
            Value::Null
        } else {
            Value::Truth(Truth::Unknown)
        });
    }
    if op == "in" || op == "not in" {
        let Value::List(list) = b else {
            bail!("in requires list");
        };
        let value = Truth::from_bool(list.contains(&a));
        return Ok(Value::Truth(if op == "in" { value } else { value.not() }));
    }
    if ["+", "-", "*", "/"].contains(&op) {
        return match (a, b) {
            (Value::Number(a), Value::Number(b)) => {
                if op == "/" && b == 0.0 {
                    return Ok(Value::Null);
                }
                let value = match op {
                    "+" => a + b,
                    "-" => a - b,
                    "*" => a * b,
                    _ => a / b,
                };
                if !value.is_finite() {
                    bail!("non-finite arithmetic result");
                }
                Ok(Value::Number(value))
            }
            (Value::Count(a), Value::Count(b)) => {
                let value = match op {
                    "+" => a.checked_add(b),
                    "-" => a.checked_sub(b),
                    "*" => a.checked_mul(b),
                    _ => a.checked_div(b),
                };
                Ok(Value::Count(value.context(
                    "count arithmetic overflow or division by zero",
                )?))
            }
            _ => bail!("arithmetic requires matching numeric types"),
        };
    }
    let ord = match (&a, &b) {
        (Value::Number(a), Value::Number(b)) => a.partial_cmp(b),
        (Value::Count(a), Value::Count(b)) => Some(a.cmp(b)),
        (Value::String(a), Value::String(b)) => Some(a.cmp(b)),
        (Value::Truth(Truth::Unknown), _) | (_, Value::Truth(Truth::Unknown)) => {
            return Ok(Value::Truth(Truth::Unknown))
        }
        (Value::Truth(a), Value::Truth(b)) => Some(a.cmp(b)),
        _ => None,
    }
    .context("comparison requires matching scalar types")?;
    Ok(Value::Truth(Truth::from_bool(match op {
        "==" => ord.is_eq(),
        "!=" => !ord.is_eq(),
        "<" => ord.is_lt(),
        "<=" => !ord.is_gt(),
        ">" => ord.is_gt(),
        ">=" => !ord.is_lt(),
        _ => bail!("unknown operator {op}"),
    })))
}
