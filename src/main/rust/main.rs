
use clap::{Parser, Subcommand};
use log::debug;

mod file;
mod helper;
mod ind;

#[derive(Subcommand)]
enum Commands {
    Create {
        #[clap(value_parser)]
        index_length: u64,
        // We explicitly use std::primitive::bool to prevent clap
        // treating this as a flag. Now it requires a value you
        // can supply false
        #[arg(long, default_value_t = true)]
        skip_first: std::primitive::bool
    },
    Filter {
        #[arg(long, default_value_t = 2)]
        cutoff: usize,
        #[clap(value_delimiter = ' ', num_args = 1..)]
        index_file_paths: Vec<String>
    }
}

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {

    #[command(subcommand)]
    command: Commands, 
}


fn main () {
    env_logger::init ();

    debug! ("hello");
    let args = Cli::parse ();

    let nucleotides = vec![
                String::from ("A"),
                String::from ("C"),
                String::from ("G"),
                String::from ("T"),
            ];

    match args.command
    {
        Commands::Create { index_length, skip_first } => {
            debug! ("create");
            ind::create_indices (nucleotides, index_length, skip_first).expect ("Failed to create indices");
        },
        Commands::Filter { cutoff, index_file_paths } => {
            debug! ("ifp: {:?}", index_file_paths);
            ind::filter_indices (nucleotides, index_file_paths, cutoff).expect ("Failed to filter indices");
        }
    }
}

