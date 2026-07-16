use clap::{Args, Parser, Subcommand};
use helicase::{Config, FastxParser, HelicaseParser, ParserOptions, input::*};
use lexichash::{LexicSketch, PartialSketch, SketchBuilder};
use std::thread;

const CONFIG: Config = ParserOptions::default()
    .ignore_headers()
    .dna_packed()
    .split_non_actg()
    .return_record(false)
    .config();

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Build a sketch from a FASTA/Q file
    Build(BuildArgs),
    /// Compare two sketches on disk
    Compare(CompareArgs),
}

#[derive(Args, Debug)]
struct BuildArgs {
    /// Input file (FASTA/Q, possibly compressed)
    input: String,
    /// Output file
    #[arg(short, long)]
    output: String,
    /// K-mer size
    #[arg(short)]
    k: usize,
    /// Prefix size
    #[arg(short, long)]
    prefix_size: usize,
    /// Number of threads [default: all]
    #[arg(short, long)]
    threads: Option<usize>,
    /// Use canonical k-mers
    #[arg(short, long)]
    canonical: bool,
}

#[derive(Args, Debug)]
struct CompareArgs {
    /// First sketch file
    sketch_1: String,
    /// Second sketch file
    sketch_2: String,
    /// Number of threads [default: all]
    #[arg(short, long)]
    threads: Option<usize>,
}

fn main() {
    let args = Cli::parse();
    match args.command {
        // Build Sketch
        Command::Build(args) => {
            let threads = args.threads.unwrap_or_else(|| {
                thread::available_parallelism()
                    .map(|n| n.get())
                    .unwrap_or(4)
            });

            // Create a parser with the desired options
            let mut parser =
                FastxParser::<CONFIG>::from_file_mmap(&args.input).expect("Cannot open the file");

            // Create sketch builder
            let builder = SketchBuilder::new(args.k, args.prefix_size, threads);
            let mut partial_sketch = PartialSketch::new(args.k, args.prefix_size);

            // Iterate over records
            while let Some(_event) = parser.next() {
                let dna = parser.get_dna_packed();
                if dna.len() >= args.k {
                    // builder.build_with(dna, &mut partial_sketch);
                    // builder.build_with_advanced::<false, false>(dna, &mut partial_sketch);
                    builder.build_with_advanced::<true, false>(dna, &mut partial_sketch);
                }
            }
            let sketch = partial_sketch.merge();
            sketch.serialize(&args.output);
        }

        // Compare two sketches
        Command::Compare(args) => {
            let sketch1 = LexicSketch::deserialize(args.sketch_1);
            let sketch2 = LexicSketch::deserialize(args.sketch_2);
            let mean = sketch1.average_match_size(&sketch2);
            let mut_rate = sketch1.get_divergence_from_mean(mean);
            println!("The average score between the two sketches is {}", mean);
            println!("Estimated mutation rate: {}%", mut_rate * 100.);
        }
    }
}
