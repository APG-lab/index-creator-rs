
use clap::{Parser, Subcommand};
use log::debug;
use std::str;

mod file;
mod helper;
mod ind;

#[derive(Subcommand)]
enum Commands {
    Create {
        #[clap(value_parser)]
        adapter_five_prime: String,
        #[clap(value_parser)]
        adapter_three_prime: String,
        #[clap(value_parser)]
        index_length: u64,
        // We explicitly use std::primitive::bool to prevent clap
        // treating this as a flag. Now it requires a value you
        // can supply false
        #[arg(long, default_value_t = true)]
        skip_first: std::primitive::bool
    },
    CreateFlorian {
        #[clap(value_parser)]
        adapter_five_prime: String,
        #[clap(value_parser)]
        adapter_three_prime: String,
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
        Commands::Create { adapter_five_prime, adapter_three_prime, index_length, skip_first } => {
            debug! ("create");
            if adapter_five_prime.bytes ().all (|x| { let xb = [x];nucleotides.contains (&str::from_utf8 (&xb).unwrap ().to_owned ()) }) &&
                adapter_three_prime.bytes ().all (|x| { let xb = [x];nucleotides.contains (&str::from_utf8 (&xb).unwrap ().to_owned ()) })
            {
                let all = ind::create_indices (&nucleotides, index_length);
                let picks = ind::pick_indices (adapter_three_prime.clone (), &nucleotides, skip_first, &all).expect ("Failed to pick indices");
                ind::output_indices (adapter_five_prime, adapter_three_prime, &nucleotides, index_length, &picks);
            }
            else
            {
                eprintln! ("Adapter five prime or adapter three prime contains chars not found in nucleotides");
            }
        },
        Commands::CreateFlorian { adapter_five_prime, adapter_three_prime, index_length, skip_first } => {
            debug! ("create");
            if adapter_five_prime.bytes ().all (|x| { let xb = [x];nucleotides.contains (&str::from_utf8 (&xb).unwrap ().to_owned ()) }) &&
                adapter_three_prime.bytes ().all (|x| { let xb = [x];nucleotides.contains (&str::from_utf8 (&xb).unwrap ().to_owned ()) })
            {
                let all = ind::create_indices (&nucleotides, index_length);
                let picks = ind::pick_indices_florian (adapter_three_prime.clone (), &nucleotides, skip_first, &all).expect ("Failed to pick indices");
                ind::output_indices (adapter_five_prime, adapter_three_prime, &nucleotides, index_length, &picks);
            }
            else
            {
                eprintln! ("Adapter five prime or adapter three prime contains chars not found in nucleotides");
            }
        },
        Commands::Filter { cutoff, index_file_paths } => {
            debug! ("ifp: {:?}", index_file_paths);
            ind::filter_indices (nucleotides, index_file_paths, cutoff).expect ("Failed to filter indices");
        }
    }
}

