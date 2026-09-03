//! Minimal SVG output for `--plot` flags: genome-track drawing with zero dependencies, styled
//! to match the paper's figures (same palette, same restrained grid). Deliberately not a chart
//! library — every plot here is genomics-shaped (coverage areas, junction arcs, site lollipops),
//! which is exactly what general chart crates lack.

use std::fmt::Write as _;

pub const BLUE: &str = "#1a6faf";
pub const DBLUE: &str = "#134f7c";
pub const LBLUE: &str = "#7fb3d5";
pub const ORANGE: &str = "#c98f1a";
pub const GRAY: &str = "#8a8a8a";
pub const LGRAY: &str = "#d9d9d9";
pub const INK: &str = "#333333";

pub struct Svg {
    pub w: f64,
    pub h: f64,
    body: String,
}

impl Svg {
    pub fn new(w: f64, h: f64) -> Svg {
        Svg { w, h, body: String::new() }
    }

    pub fn rect(&mut self, x: f64, y: f64, w: f64, h: f64, fill: &str, opacity: f64) {
        let _ = writeln!(
            self.body,
            r#"<rect x="{x:.1}" y="{y:.1}" width="{w:.1}" height="{h:.1}" fill="{fill}" opacity="{opacity}"/>"#
        );
    }

    pub fn line(&mut self, x1: f64, y1: f64, x2: f64, y2: f64, stroke: &str, width: f64) {
        let _ = writeln!(
            self.body,
            r#"<line x1="{x1:.1}" y1="{y1:.1}" x2="{x2:.1}" y2="{y2:.1}" stroke="{stroke}" stroke-width="{width}"/>"#
        );
    }

    pub fn path(&mut self, d: &str, fill: &str, stroke: &str, width: f64, opacity: f64) {
        let _ = writeln!(
            self.body,
            r#"<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{width}" opacity="{opacity}" stroke-linecap="round"/>"#
        );
    }

    pub fn circle(&mut self, cx: f64, cy: f64, r: f64, fill: &str, stroke: &str, width: f64) {
        let _ = writeln!(
            self.body,
            r#"<circle cx="{cx:.1}" cy="{cy:.1}" r="{r:.1}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>"#
        );
    }

    /// anchor: "start" | "middle" | "end"
    pub fn text(&mut self, x: f64, y: f64, s: &str, size: f64, fill: &str, anchor: &str) {
        let esc = s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;");
        let _ = writeln!(
            self.body,
            r#"<text x="{x:.1}" y="{y:.1}" font-size="{size}" fill="{fill}" text-anchor="{anchor}">{esc}</text>"#
        );
    }

    pub fn text_italic(&mut self, x: f64, y: f64, s: &str, size: f64, fill: &str, anchor: &str) {
        let esc = s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;");
        let _ = writeln!(
            self.body,
            r#"<text x="{x:.1}" y="{y:.1}" font-size="{size}" fill="{fill}" text-anchor="{anchor}" font-style="italic">{esc}</text>"#
        );
    }

    pub fn finish(self) -> String {
        format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {:.0} {:.0}\" \
             font-family=\"Helvetica,Arial,sans-serif\">\n<rect width=\"{:.0}\" height=\"{:.0}\" \
             fill=\"#ffffff\"/>\n{}</svg>\n",
            self.w, self.h, self.w, self.h, self.body
        )
    }
}

/// Round a span to a "nice" tick step (1/2/5 × 10^k) giving 4–8 ticks.
pub fn nice_step(span: f64) -> f64 {
    let raw = span / 6.0;
    let mag = 10f64.powf(raw.log10().floor());
    let n = raw / mag;
    let m = if n < 1.5 {
        1.0
    } else if n < 3.5 {
        2.0
    } else if n < 7.5 {
        5.0
    } else {
        10.0
    };
    m * mag
}

/// Genomic x-axis with ticks in kb/Mb units.
pub fn genomic_axis(svg: &mut Svg, x0: f64, x1: f64, y: f64, gstart: f64, gend: f64, chrom: &str) {
    svg.line(x0, y, x1, y, "#777777", 0.9);
    let step = nice_step(gend - gstart);
    let sx = |g: f64| x0 + (g - gstart) / (gend - gstart) * (x1 - x0);
    let mut t = (gstart / step).ceil() * step;
    while t <= gend {
        let x = sx(t);
        svg.line(x, y, x, y + 4.0, "#777777", 0.9);
        let label = if gend >= 1e6 {
            format!("{:.2} Mb", t / 1e6)
        } else {
            format!("{:.0} kb", t / 1e3)
        };
        svg.text(x, y + 15.0, &label, 9.0, "#555555", "middle");
        t += step;
    }
    svg.text((x0 + x1) / 2.0, y + 29.0, chrom, 10.0, INK, "middle");
}
