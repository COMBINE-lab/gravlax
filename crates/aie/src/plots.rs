//! The `--plot` renderers: sashimi-style region portraits, APA site lollipops, and the EM
//! reliability diagram. Everything draws through `viz` (plain SVG, paper palette); everything
//! takes data the command already computed.

use crate::viz::{self, Svg};
use anyhow::{Context, Result};
use std::path::Path;

/// Write a finished SVG string to `out`; a `.png` extension rasterizes it via resvg at 2x.
fn write_plot(out: &Path, svg: String) -> Result<()> {
    if out.extension().is_some_and(|e| e.eq_ignore_ascii_case("png")) {
        let mut opt = resvg::usvg::Options::default();
        {
            let db = opt.fontdb_mut();
            db.load_system_fonts();
            // Headless hosts often have no fonts at all; ship one (DejaVu Sans, free license).
            db.load_font_data(include_bytes!("../assets/DejaVuSans.ttf").to_vec());
            db.set_sans_serif_family("DejaVu Sans");
        }
        opt.font_family = "DejaVu Sans".to_string();
        let tree = resvg::usvg::Tree::from_str(&svg, &opt).context("parsing generated SVG")?;
        let size = tree.size();
        let (w, h) = ((size.width() * 2.0) as u32, (size.height() * 2.0) as u32);
        let mut pixmap = resvg::tiny_skia::Pixmap::new(w, h).context("pixmap")?;
        resvg::render(
            &tree,
            resvg::tiny_skia::Transform::from_scale(2.0, 2.0),
            &mut pixmap.as_mut(),
        );
        pixmap.save_png(out)?;
    } else {
        std::fs::write(out, svg)?;
    }
    Ok(())
}

pub struct GeneBox {
    pub name: String,
    pub start: f64,
    pub end: f64,
    pub rev: bool,
}

/// Sashimi-style portrait of a window: per-strand molecule coverage with junction arcs whose
/// weight tracks molecule support, plus an optional gene underlay.
#[allow(clippy::too_many_arguments)]
pub fn region_plot(
    out: &Path,
    chrom: &str,
    start: u32,
    end: u32,
    cov: [&[u32]; 2],                      // [+, −] per-base molecule coverage
    juncs: [&[(u32, u32, u64)]; 2],        // [+, −] (donor, acceptor, molecules)
    genes: &[GeneBox],
    title: &str,
) -> Result<()> {
    let (w, x0, x1) = (900.0, 60.0, 870.0);
    let span = (end - start) as f64;
    let sx = |g: f64| x0 + (g - start as f64) / span * (x1 - x0);
    let both = !cov[0].is_empty() && cov[0].iter().any(|&v| v > 0) && cov[1].iter().any(|&v| v > 0);
    let panel_h = 130.0;
    let n_panels = if both { 2 } else { 1 };
    let gene_h = if genes.is_empty() { 0.0 } else { 34.0 };
    let h = 34.0 + n_panels as f64 * panel_h + gene_h + 44.0;
    let mut svg = Svg::new(w, h);
    svg.text(x0, 18.0, title, 11.5, viz::INK, "start");

    let mut y_top = 30.0;
    for (si, (cv, js)) in cov.iter().zip(juncs.iter()).enumerate() {
        if cv.iter().all(|&v| v == 0) {
            continue;
        }
        let color = if si == 0 { viz::BLUE } else { viz::ORANGE };
        let base = y_top + panel_h - 26.0;
        let max = cv.iter().copied().max().unwrap_or(1).max(1) as f64;
        // Coverage area, downsampled to pixel columns.
        let mut d = format!("M{:.1},{:.1}", x0, base);
        let px = (x1 - x0) as usize;
        for p in 0..=px {
            let lo = ((p as f64 / px as f64 * cv.len() as f64) as usize).min(cv.len() - 1);
            let hi = ((((p + 1) as f64 / px as f64 * cv.len() as f64) as usize).max(lo + 1)).min(cv.len());
            let v = cv[lo..hi].iter().copied().max().unwrap_or(0) as f64;
            let y = base - (v / max) * (panel_h - 52.0);
            d.push_str(&format!(" L{:.1},{:.1}", x0 + p as f64, y));
        }
        d.push_str(&format!(" L{x1:.1},{base:.1} Z"));
        svg.path(&d, color, "none", 0.0, 0.30);
        svg.line(x0, base, x1, base, "#999999", 0.8);
        svg.text(x0 - 6.0, y_top + 12.0, &format!("{} strand", if si == 0 { "+" } else { "−" }), 9.5, viz::INK, "start");
        svg.text(x0 - 6.0, base + 3.0, "0", 8.0, "#777777", "end");
        svg.text(x0 - 6.0, base - (panel_h - 52.0) + 3.0, &format!("{}", max as u64), 8.0, "#777777", "end");
        // Junction arcs, heaviest on top of the draw order; label the strongest few.
        let mut js2: Vec<&(u32, u32, u64)> = js.iter().filter(|j| j.2 >= 2).collect();
        js2.sort_unstable_by_key(|j| j.2);
        if js2.len() > 80 {
            let cut = js2.len() - 80;
            js2.drain(..cut);
        }
        let jmax = js2.last().map(|j| j.2).unwrap_or(1).max(1) as f64;
        let n_label = 6.min(js2.len());
        for (rank, (dn, ac, cnt)) in js2.iter().enumerate() {
            let (xa, xb) = (sx(*dn as f64), sx(*ac as f64));
            if xb - xa < 2.0 {
                continue;
            }
            let frac = (*cnt as f64).ln_1p() / jmax.ln_1p();
            let peak = base - (panel_h - 66.0) * (0.35 + 0.65 * frac) - 8.0;
            let width = 0.8 + 2.6 * frac;
            svg.path(
                &format!("M{xa:.1},{base:.1} Q{:.1},{peak:.1} {xb:.1},{base:.1}", (xa + xb) / 2.0),
                "none",
                if si == 0 { viz::DBLUE } else { "#8a6110" },
                width,
                0.85,
            );
            if rank + n_label >= js2.len() {
                svg.text((xa + xb) / 2.0, peak - 3.0, &format!("{cnt}"), 8.0, viz::INK, "middle");
            }
        }
        y_top += panel_h;
    }
    // Gene underlay.
    if !genes.is_empty() {
        let gy = y_top + 8.0;
        for g in genes {
            let (a, b) = (sx(g.start.max(start as f64)), sx(g.end.min(end as f64)));
            svg.rect(a, gy, (b - a).max(2.0), 7.0, viz::GRAY, 0.8);
            let label = format!("{}{}", g.name, if g.rev { " ◂" } else { " ▸" });
            svg.text(((a + b) / 2.0).clamp(x0 + 15.0, x1 - 15.0), gy + 18.0, &label, 8.5, "#555555", "middle");
        }
        y_top += gene_h;
    }
    viz::genomic_axis(&mut svg, x0, x1, y_top + 6.0, start as f64, end as f64, chrom);
    write_plot(out, svg.finish())?;
    Ok(())
}

pub struct SiteDot {
    pub cp: u32,
    pub rev: bool,
    pub umis: usize,
    pub ip: bool,
    pub gc: Vec<u64>,
}

/// 3'-site lollipops on a genomic axis; flagged internal-priming sites drawn hollow gray. With
/// groups, a per-site usage-proportion band sits under the axis.
pub fn apa_plot(
    out: &Path,
    chrom: &str,
    start: u32,
    end: u32,
    sites: &[SiteDot],
    group_names: &[String],
    title: &str,
) -> Result<()> {
    let (w, x0, x1) = (900.0, 60.0, 870.0);
    let span = (end - start) as f64;
    let sx = |g: f64| x0 + (g - start as f64) / span * (x1 - x0);
    let band_h = if group_names.len() >= 2 { 42.0 } else { 0.0 };
    let plot_h = 150.0;
    let h = 34.0 + plot_h + 30.0 + band_h + 42.0;
    let mut svg = Svg::new(w, h);
    svg.text(x0, 18.0, title, 11.5, viz::INK, "start");
    let base = 30.0 + plot_h;
    let umax = sites.iter().map(|s| s.umis).max().unwrap_or(1).max(1) as f64;
    svg.line(x0, base, x1, base, "#999999", 0.8);
    let sy = |u: usize| base - (u as f64).ln_1p() / umax.ln_1p() * (plot_h - 18.0);
    // Reference stems at powers of ten on the left.
    for refv in [10usize, 100, 1000, 10000] {
        if (refv as f64) < umax * 1.2 {
            svg.text(x0 - 6.0, sy(refv) + 3.0, &format!("{refv}"), 8.0, "#777777", "end");
            svg.line(x0 - 2.0, sy(refv), x0 + 2.0, sy(refv), "#777777", 0.8);
        }
    }
    svg.text(x0 - 34.0, 30.0 + plot_h / 2.0, "UMIs", 9.0, "#555555", "middle");
    for s in sites {
        let x = sx(s.cp as f64);
        let y = sy(s.umis);
        if s.ip {
            svg.line(x, base, x, y, viz::LGRAY, 1.0);
            svg.circle(x, y, 3.2, "#ffffff", viz::GRAY, 1.1);
        } else {
            let c = if s.rev { viz::ORANGE } else { viz::BLUE };
            svg.line(x, base, x, y, c, 1.2);
            svg.circle(x, y, 3.6, c, "none", 0.0);
        }
    }
    // Legend.
    svg.circle(x1 - 190.0, 16.0, 3.4, viz::BLUE, "none", 0.0);
    svg.text(x1 - 183.0, 19.0, "+ site", 8.5, viz::INK, "start");
    svg.circle(x1 - 140.0, 16.0, 3.4, viz::ORANGE, "none", 0.0);
    svg.text(x1 - 133.0, 19.0, "− site", 8.5, viz::INK, "start");
    svg.circle(x1 - 88.0, 16.0, 3.2, "#ffffff", viz::GRAY, 1.1);
    svg.text(x1 - 81.0, 19.0, "internal priming", 8.5, viz::INK, "start");
    // Group-proportion band.
    let mut y_axis = base + 8.0;
    if band_h > 0.0 {
        let by = base + 12.0;
        let palette = [viz::BLUE, viz::ORANGE, "#5a9e6f", "#9467bd", "#8c564b", "#7f7f7f"];
        for s in sites.iter().filter(|s| !s.ip) {
            let tot: u64 = s.gc.iter().sum();
            if tot < 10 {
                continue;
            }
            let x = sx(s.cp as f64);
            let mut y = by;
            for (gi, &n) in s.gc.iter().enumerate() {
                let hh = n as f64 / tot as f64 * (band_h - 14.0);
                svg.rect(x - 2.2, y, 4.4, hh, palette[gi % palette.len()], 0.9);
                y += hh;
            }
        }
        for (gi, gname) in group_names.iter().enumerate() {
            let lx = x0 + 90.0 * gi as f64;
            svg.rect(lx, by + band_h - 10.0, 8.0, 8.0, palette[gi % palette.len()], 0.9);
            svg.text(lx + 12.0, by + band_h - 3.0, gname, 8.5, viz::INK, "start");
        }
        svg.text_italic(x0, by - 2.0, "usage share by population", 8.0, "#777777", "start");
        y_axis = by + band_h + 4.0;
    }
    viz::genomic_axis(&mut svg, x0, x1, y_axis, start as f64, end as f64, chrom);
    write_plot(out, svg.finish())?;
    Ok(())
}

/// Reliability diagram for the EM masked-recovery run: per mode, empirical accuracy in each
/// max-responsibility decile against the diagonal ("r means what it says").
pub fn em_plot(
    out: &Path,
    cals: &[(String, [[u64; 2]; 10], f64, u64)],
    title: &str,
) -> Result<()> {
    let (w, h, x0, y0, x1, y1) = (420.0, 400.0, 60.0, 40.0, 390.0, 330.0);
    let mut svg = Svg::new(w, h);
    svg.text(x0, 20.0, title, 11.5, viz::INK, "start");
    for f in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let x = x0 + f * (x1 - x0);
        let y = y1 - f * (y1 - y0);
        svg.line(x, y1, x, y1 + 4.0, "#777777", 0.8);
        svg.text(x, y1 + 15.0, &format!("{f:.2}"), 8.5, "#555555", "middle");
        svg.line(x0 - 4.0, y, x0, y, "#777777", 0.8);
        svg.text(x0 - 7.0, y + 3.0, &format!("{f:.2}"), 8.5, "#555555", "end");
        svg.line(x0, y, x1, y, "#eeeeee", 0.7);
    }
    svg.path(&format!("M{x0},{y1} L{x1},{y0}"), "none", "#bbbbbb", 1.0, 1.0);
    svg.line(x0, y1, x1, y1, "#777777", 0.9);
    svg.line(x0, y0, x0, y1, "#777777", 0.9);
    svg.text((x0 + x1) / 2.0, y1 + 30.0, "max responsibility (predicted probability)", 9.5, viz::INK, "middle");
    svg.text_italic(18.0, (y0 + y1) / 2.0, "empirical accuracy", 9.5, viz::INK, "middle");
    let palette = [viz::GRAY, viz::LBLUE, viz::BLUE];
    for (mi, (name, cal, top1, n)) in cals.iter().enumerate() {
        let color = palette[mi % palette.len()];
        let mut prev: Option<(f64, f64)> = None;
        for (b, [cn, cc]) in cal.iter().enumerate() {
            if *cn == 0 {
                continue;
            }
            let x = x0 + (b as f64 + 0.5) / 10.0 * (x1 - x0);
            let y = y1 - (*cc as f64 / *cn as f64) * (y1 - y0);
            let r = 2.0 + (*cn as f64).ln_1p() * 0.55;
            if let Some((px, py)) = prev {
                svg.line(px, py, x, y, color, 1.2);
            }
            svg.circle(x, y, r, color, "#ffffff", 0.7);
            prev = Some((x, y));
        }
        let ly = y0 + 14.0 * mi as f64;
        svg.circle(x0 + 12.0, ly, 4.0, color, "none", 0.0);
        svg.text(x0 + 20.0, ly + 3.0, &format!("{name} (top-1 {:.1}%, n={n})", 100.0 * top1), 8.5, viz::INK, "start");
    }
    write_plot(out, svg.finish())?;
    Ok(())
}
