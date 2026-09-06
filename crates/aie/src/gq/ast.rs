//! Syntax is resource-independent. Byte spans survive binding and type checking.
use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct Expr {
    pub at: usize,
    pub kind: Kind,
}

#[derive(Clone, Debug, Serialize)]
pub enum Kind {
    Name(String),
    Resource(String),
    String(String),
    Number(f64),
    Count(u64),
    Bases(u32),
    Bool(bool),
    Null,
    Unknown,
    Locus {
        junction: bool,
        chrom: String,
        start: u32,
        end: u32,
        reverse: Option<bool>,
    },
    List(Vec<Expr>),
    Unary(String, Box<Expr>),
    Binary(String, Box<Expr>, Box<Expr>),
    Call(String, Vec<(Option<String>, Expr)>),
    Field(Box<Expr>, String),
    Quant {
        all: bool,
        domain: String,
        body: Box<Expr>,
    },
    Path {
        order: String,
        body: Box<Expr>,
    },
}

#[derive(Clone, Debug, Serialize)]
pub struct Function {
    pub name: String,
    pub parameters: Vec<(String, Option<String>)>,
    pub returns: Option<String>,
    pub exported: bool,
    pub body: Expr,
}

#[derive(Clone, Debug, Serialize)]
pub enum Stage {
    Support(String, Vec<(String, Expr)>),
    Within(Expr),
    Where(Expr),
    Derive(Vec<(String, Expr)>),
    Tally(Vec<(String, Expr)>, Vec<(String, Expr)>),
    Summarize(Vec<(String, Expr)>, Vec<(String, Expr)>),
    Select(Vec<(String, Expr)>),
    Sort(Vec<(String, Expr)>),
    Take(usize),
}

#[derive(Clone, Debug, Serialize)]
pub struct Document {
    pub header: Vec<(String, Expr)>,
    pub bindings: Vec<(String, Expr)>,
    pub functions: Vec<Function>,
    pub source: Expr,
    pub stages: Vec<Stage>,
}
