//! Expression parser
//!
//! Converts string expressions to Expr trees

use super::expr::Expr;
use super::{SymbolicError, SymbolicResult};

/// Parse a string into an expression
pub fn parse(input: &str) -> SymbolicResult<Expr> {
    let tokens = tokenize(input)?;
    let mut parser = Parser::new(tokens);
    parser.parse_expression()
}

#[derive(Debug, Clone, PartialEq)]
enum Token {
    Number(i64),
    Symbol(String),
    Plus,
    Minus,
    Star,
    Slash,
    Caret,
    LParen,
    RParen,
    Comma,
}

struct Parser {
    tokens: Vec<Token>,
    pos: usize,
}

impl Parser {
    fn new(tokens: Vec<Token>) -> Self {
        Parser { tokens, pos: 0 }
    }

    fn current(&self) -> Option<&Token> {
        self.tokens.get(self.pos)
    }

    fn advance(&mut self) {
        if self.pos < self.tokens.len() {
            self.pos += 1;
        }
    }

    fn parse_expression(&mut self) -> SymbolicResult<Expr> {
        self.parse_additive()
    }

    fn parse_additive(&mut self) -> SymbolicResult<Expr> {
        let mut left = self.parse_multiplicative()?;

        while let Some(token) = self.current() {
            match token {
                Token::Plus => {
                    self.advance();
                    let right = self.parse_multiplicative()?;
                    left = Expr::add(left, right);
                }
                Token::Minus => {
                    self.advance();
                    let right = self.parse_multiplicative()?;
                    left = Expr::add(left, Expr::mul(Expr::num(-1), right));
                }
                _ => break,
            }
        }

        Ok(left)
    }

    fn parse_multiplicative(&mut self) -> SymbolicResult<Expr> {
        let mut left = self.parse_power()?;

        while let Some(token) = self.current() {
            match token {
                Token::Star => {
                    self.advance();
                    let right = self.parse_power()?;
                    left = Expr::mul(left, right);
                }
                Token::Slash => {
                    self.advance();
                    let right = self.parse_power()?;
                    left = Expr::mul(left, Expr::pow(right, Expr::num(-1)));
                }
                _ => break,
            }
        }

        Ok(left)
    }

    fn parse_power(&mut self) -> SymbolicResult<Expr> {
        let mut left = self.parse_primary()?;

        if let Some(Token::Caret) = self.current() {
            self.advance();
            let right = self.parse_power()?;
            left = Expr::pow(left, right);
        }

        Ok(left)
    }

    fn parse_primary(&mut self) -> SymbolicResult<Expr> {
        match self.current() {
            Some(Token::Number(n)) => {
                let n = *n;
                self.advance();
                Ok(Expr::num(n))
            }
            Some(Token::Symbol(s)) => {
                let s = s.clone();
                self.advance();

                // Check if this is a function call
                if let Some(Token::LParen) = self.current() {
                    self.advance();
                    let mut args = Vec::new();

                    if !matches!(self.current(), Some(Token::RParen)) {
                        args.push(self.parse_expression()?);

                        while let Some(Token::Comma) = self.current() {
                            self.advance();
                            args.push(self.parse_expression()?);
                        }
                    }

                    if !matches!(self.current(), Some(Token::RParen)) {
                        return Err(SymbolicError::ParseError("Expected ')'".to_string()));
                    }
                    self.advance();

                    Ok(Expr::func(s, args))
                } else {
                    Ok(Expr::sym(s))
                }
            }
            Some(Token::LParen) => {
                self.advance();
                let expr = self.parse_expression()?;
                if !matches!(self.current(), Some(Token::RParen)) {
                    return Err(SymbolicError::ParseError("Expected ')'".to_string()));
                }
                self.advance();
                Ok(expr)
            }
            Some(Token::Minus) => {
                self.advance();
                let expr = self.parse_primary()?;
                Ok(Expr::mul(Expr::num(-1), expr))
            }
            _ => Err(SymbolicError::ParseError(format!(
                "Unexpected token: {:?}",
                self.current()
            ))),
        }
    }
}

fn tokenize(input: &str) -> SymbolicResult<Vec<Token>> {
    let mut tokens = Vec::new();
    let chars: Vec<char> = input.chars().collect();
    let mut i = 0;

    while i < chars.len() {
        match chars[i] {
            ' ' | '\t' | '\n' => {
                i += 1;
            }
            '+' => {
                tokens.push(Token::Plus);
                i += 1;
            }
            '-' => {
                tokens.push(Token::Minus);
                i += 1;
            }
            '*' => {
                tokens.push(Token::Star);
                i += 1;
            }
            '/' => {
                tokens.push(Token::Slash);
                i += 1;
            }
            '^' => {
                tokens.push(Token::Caret);
                i += 1;
            }
            '(' => {
                tokens.push(Token::LParen);
                i += 1;
            }
            ')' => {
                tokens.push(Token::RParen);
                i += 1;
            }
            ',' => {
                tokens.push(Token::Comma);
                i += 1;
            }
            c if c.is_ascii_digit() => {
                let start = i;
                while i < chars.len() && chars[i].is_ascii_digit() {
                    i += 1;
                }
                let num_str: String = chars[start..i].iter().collect();
                let num = num_str
                    .parse::<i64>()
                    .map_err(|e| SymbolicError::ParseError(format!("Invalid number: {}", e)))?;
                tokens.push(Token::Number(num));
            }
            c if c.is_alphabetic() => {
                let start = i;
                while i < chars.len() && (chars[i].is_alphanumeric() || chars[i] == '_') {
                    i += 1;
                }
                let symbol: String = chars[start..i].iter().collect();
                tokens.push(Token::Symbol(symbol));
            }
            c => {
                return Err(SymbolicError::ParseError(format!(
                    "Unexpected character: {}",
                    c
                )));
            }
        }
    }

    Ok(tokens)
}

