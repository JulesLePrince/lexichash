// Benchmarks build throughput and compare time across prefix sizes.
// Both input files are loaded into RAM upfront
// so the timed sections measure sketch building/comparison, not disk reads.
//
// Run e.g. `cargo run --release --example benchmark -- chm13v2.0.fa chr1.fa`.
use clap::Parser;
use core::hint::black_box;
use helicase::{Config, FastxParser, HelicaseParser, ParserOptions, input::*};
use lexichash::{LexicSketch, PartialSketch, SketchBuilder};
use std::time::Instant;

const CONFIG: Config = ParserOptions::default()
    .ignore_headers()
    .dna_packed()
    .split_non_actg()
    .return_record(false)
    .config();

#[derive(Parser, Debug)]
#[command(about = "Benchmark lexichash build/compare across prefix sizes")]
struct Cli {
    /// First FASTA/Q file (build throughput benchmark, one side of the compare pair)
    fasta_1: String,
    /// Second FASTA/Q file (other side of the compare pair)
    fasta_2: String,
    /// K-mer size
    #[arg(short, default_value_t = 19)]
    k: usize,
    /// Minimum prefix size (inclusive)
    #[arg(long, default_value_t = 4)]
    p_min: usize,
    /// Maximum prefix size (inclusive)
    #[arg(long, default_value_t = 8)]
    p_max: usize,
    /// Number of threads
    #[arg(short, long, default_value_t = 1)]
    threads: usize,
    /// Timed runs for the build benchmark
    #[arg(long, default_value_t = 3)]
    build_runs: usize,
    /// Warmup runs for the build benchmark
    #[arg(long, default_value_t = 1)]
    build_warmup: usize,
    /// Timed runs for the compare benchmark
    #[arg(long, default_value_t = 10)]
    compare_runs: usize,
    /// Warmup runs for the compare benchmark
    #[arg(long, default_value_t = 2)]
    compare_warmup: usize,
    /// Emit a Typst table instead of plain text
    #[arg(long)]
    typst: bool,
}

fn build_sketch(data: &[u8], k: usize, prefix_size: usize, threads: usize) -> LexicSketch {
    let mut parser = FastxParser::<CONFIG>::from_slice(data).expect("Cannot parse the data");
    let builder = SketchBuilder::new(k, prefix_size, threads);
    let mut partial_sketch = PartialSketch::new(k, prefix_size);
    while let Some(_event) = parser.next() {
        let dna = parser.get_dna_packed();
        if dna.len() >= k {
            builder.process_seq(dna, &mut partial_sketch);
        }
    }
    partial_sketch.merge()
}

/// Total number of bases actually fed to the sketch builder (matches the
/// `dna.len() >= k` filter in `build_sketch`), used for the Gbp/s throughput.
fn count_bases(data: &[u8], k: usize) -> usize {
    let mut parser = FastxParser::<CONFIG>::from_slice(data).expect("Cannot parse the data");
    let mut bases = 0;
    while let Some(_event) = parser.next() {
        let dna = parser.get_dna_packed();
        if dna.len() >= k {
            bases += dna.len();
        }
    }
    bases
}

/// Formats `value` with `sig_figs` significant digits (leading zeros after
/// the decimal point don't count).
fn format_sig_figs(value: f64, sig_figs: usize) -> String {
    if value == 0.0 || !value.is_finite() {
        return format!("{value}");
    }
    let magnitude = value.abs().log10().floor() as i32;
    let decimals = sig_figs as i32 - 1 - magnitude;
    let factor = 10f64.powi(decimals);
    let rounded = (value * factor).round() / factor;
    let display_decimals = decimals.max(0) as usize;
    format!("{rounded:.display_decimals$}")
}

fn mean_std(values: &[f64]) -> (f64, f64) {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
    (mean, variance.sqrt())
}

fn bench_build_throughput(
    data: &[u8],
    k: usize,
    p: usize,
    threads: usize,
    runs: usize,
    warmup: usize,
    size_gbp: f64,
) -> (f64, f64) {
    for _ in 0..warmup {
        black_box(data);
        build_sketch(data, k, p, threads);
    }
    let throughputs: Vec<f64> = (0..runs)
        .map(|_| {
            let start = Instant::now();
            let sketch = build_sketch(data, k, p, threads);
            let elapsed = start.elapsed().as_secs_f64();
            black_box(sketch);
            size_gbp / elapsed
        })
        .collect();
    mean_std(&throughputs)
}

fn bench_compare_time(
    sketch_1: &LexicSketch,
    sketch_2: &LexicSketch,
    runs: usize,
    warmup: usize,
) -> (f64, f64) {
    for _ in 0..warmup {
        black_box(sketch_1);
        black_box(sketch_2);
        sketch_1.average_match_size(sketch_2);
    }
    let times_us: Vec<f64> = (0..runs)
        .map(|_| {
            let start = Instant::now();
            let score = sketch_1.average_match_size(sketch_2);
            let elapsed = start.elapsed().as_secs_f64() * 1_000_000.0;
            black_box(score);
            elapsed
        })
        .collect();
    mean_std(&times_us)
}

fn print_text_table(p_values: &[usize], build: &[(f64, f64)], compare: &[(f64, f64)]) {
    let mut header = vec!["metric".to_string()];
    header.extend(p_values.iter().map(|p| format!("p={p}")));

    let build_row: Vec<String> = std::iter::once("build throughput (Gbp/s)".to_string())
        .chain(build.iter().map(|(mean, std)| {
            format!(
                "{} +/- {}",
                format_sig_figs(*mean, 3),
                format_sig_figs(*std, 3)
            )
        }))
        .collect();
    let compare_row: Vec<String> = std::iter::once("compare time (us)".to_string())
        .chain(compare.iter().map(|(mean, std)| {
            format!(
                "{} +/- {}",
                format_sig_figs(*mean, 3),
                format_sig_figs(*std, 3)
            )
        }))
        .collect();

    let rows = [header, build_row, compare_row];
    let cols = rows[0].len();
    let widths: Vec<usize> = (0..cols)
        .map(|c| rows.iter().map(|r| r[c].len()).max().unwrap_or(0))
        .collect();

    for row in &rows {
        let line: Vec<String> = row
            .iter()
            .zip(&widths)
            .map(|(cell, width)| format!("{cell:<width$}"))
            .collect();
        println!("{}", line.join("  "));
    }
}

fn print_typst_table(p_values: &[usize], build: &[(f64, f64)], compare: &[(f64, f64)]) {
    let cols = p_values.len() + 1;
    println!("#table(");
    println!("  columns: {cols},");
    let align = std::iter::once("left")
        .chain(std::iter::repeat_n("center", p_values.len()))
        .collect::<Vec<_>>()
        .join(", ");
    println!("  align: ({align}),");

    let header: Vec<String> = std::iter::once("[]".to_string())
        .chain(p_values.iter().map(|p| format!("[*{p}*]")))
        .collect();
    println!("  table.header{},", header.join(""));

    let build_cells: Vec<String> = std::iter::once("[*build throughput (Gbp/s)*]".to_string())
        .chain(build.iter().map(|(mean, std)| {
            format!(
                "[{} $plus.minus$ {}]",
                format_sig_figs(*mean, 3),
                format_sig_figs(*std, 3)
            )
        }))
        .collect();
    println!("  {},", build_cells.join(", "));

    let compare_cells: Vec<String> = std::iter::once("[*compare time (µs)*]".to_string())
        .chain(compare.iter().map(|(mean, std)| {
            format!(
                "[{} $plus.minus$ {}]",
                format_sig_figs(*mean, 3),
                format_sig_figs(*std, 3)
            )
        }))
        .collect();
    println!("  {},", compare_cells.join(", "));

    println!(")");
}

fn main() {
    let cli = Cli::parse();

    // Load both files into RAM upfront so the timed sections only measure
    // sketch building/comparison, not disk reads.
    let fasta_1_data = std::fs::read(&cli.fasta_1).expect("Cannot read fasta_1");
    let fasta_2_data = std::fs::read(&cli.fasta_2).expect("Cannot read fasta_2");
    let fasta_1_size_gbp = count_bases(&fasta_1_data, cli.k) as f64 / 1e9;

    let p_values: Vec<usize> = (cli.p_min..=cli.p_max).collect();

    let mut build_throughput = Vec::with_capacity(p_values.len());
    let mut compare_time = Vec::with_capacity(p_values.len());

    for &p in &p_values {
        build_throughput.push(bench_build_throughput(
            &fasta_1_data,
            cli.k,
            p,
            cli.threads,
            cli.build_runs,
            cli.build_warmup,
            fasta_1_size_gbp,
        ));

        let sketch_1 = build_sketch(&fasta_1_data, cli.k, p, cli.threads);
        let sketch_2 = build_sketch(&fasta_2_data, cli.k, p, cli.threads);
        compare_time.push(bench_compare_time(
            &sketch_1,
            &sketch_2,
            cli.compare_runs,
            cli.compare_warmup,
        ));
    }

    if cli.typst {
        print_typst_table(&p_values, &build_throughput, &compare_time);
    } else {
        print_text_table(&p_values, &build_throughput, &compare_time);
    }
}
