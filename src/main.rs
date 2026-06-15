use clap::{Args, Parser, Subcommand};
use lexichash::sketch_builder::SketchBuilder;
use std::thread;

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
        Command::Build(args) => {
            let threads = args.threads.unwrap_or_else(|| {
                thread::available_parallelism()
                    .map(|n| n.get())
                    .unwrap_or(4)
            });
            let builder = SketchBuilder::new(args.k, args.prefix_size, threads);
            let seq = todo!();
            let sketch = builder.build(&seq);
        }
        Command::Compare(args) => todo!(),
    }
}
