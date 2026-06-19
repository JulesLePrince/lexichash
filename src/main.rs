use clap::{Args, Parser, Subcommand};
use epserde::ser::Serialize;
use helicase::{Config, FastxParser, HelicaseParser, ParserOptions, input::FromFile};
use lexichash::{LexicSketch, SketchBuilder};
use std::thread;

const CONFIG: Config = ParserOptions::default().and_dna_packed().config();

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

            // Create sketch builder
            let builder = SketchBuilder::new(args.k, args.prefix_size, threads);

            // Create a parser with the desired options
            let mut parser =
                FastxParser::<CONFIG>::from_file(&args.input).expect("Cannot open the file");

            // Iterate over records
            if let Some(_event) = parser.next() {
                // get a reference to the packed sequence
                let seq = parser.get_dna_packed();
                // Build the sketch and serialize it
                let sketch = builder.build(&seq);
                sketch.serialize(args.output);
            }
        }

        // Compare two sketches
        Command::Compare(args) => {
            let sketch1 = LexicSketch::deserialize(args.sketch_1);
            let sketch2 = LexicSketch::deserialize(args.sketch_2);
            let score = sketch1.get_score(&sketch2);
            println!("The score between sketch 1 and sketch 2 is {}", score);
        },
    }
}
