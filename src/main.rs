// sablast: Spectral Burrows-Wheeler transform accelerated local alignment search
//
// Copyright 2024 Tommi Mäklin [tommi@maklin.fi].

// Copyrights in this project are retained by contributors. No copyright assignment
// is required to contribute to this project.

// Except as otherwise noted (below and/or in individual files), this
// project is licensed under the Apache License, Version 2.0
// <LICENSE-APACHE> or <http://www.apache.org/licenses/LICENSE-2.0> or
// the MIT license, <LICENSE-MIT> or <http://opensource.org/licenses/MIT>,
// at your option.
//
use std::ffi::OsString;

use clap::Parser;
use log::{info, Record, Level, Metadata};
use sbwt::SbwtIndexVariant;

// Command-line interface
mod cli;

// Subcommand implementations
mod build;
mod map;

// Logger implementation
struct Logger;

impl log::Log for Logger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= Level::Info
    }

    fn log(&self, record: &Record) {
        if self.enabled(record.metadata()) {
            println!("{} - {}", record.level(), record.args());
        }
    }

    fn flush(&self) {}
}

fn init_log(log_max_level: usize) {
    stderrlog::new()
	.module(module_path!())
	.quiet(false)
	.verbosity(log_max_level)
	.timestamp(stderrlog::Timestamp::Off)
	.init()
	.unwrap();
}

// Use `sablast` to list the available commands or `sablast <command>` to run.
fn main() {
    let cli = cli::Cli::parse();

    // Subcommands:
    match &cli.command {
        // Run the full pipeline
        Some(cli::Commands::Build {
	    seq_files,
	    input_list,
	    output_prefix,
	    num_threads,
	    mem_gb,
	    temp_dir,
	    verbose,
        }) => {
	    init_log(if *verbose { 2 } else { 1 });
	    info!("Building SBWT index...");

            let sbwt_params = build::SBWTParams {
		num_threads: *num_threads,
		mem_gb: *mem_gb,
		temp_dir: Some(std::path::PathBuf::from(OsString::from(temp_dir.clone().unwrap()))),
		index_prefix: output_prefix.clone(),
                ..Default::default()
            };

	    // TODO handle multiple files and `input_list`
	    info!("Serializing SBWT index...");
	    let (sbwt, lcs) = build::build_sbwt(&seq_files[0], &Some(sbwt_params.clone()));
	    build::serialize_sbwt(sbwt, &lcs, &Some(sbwt_params));

	},
        Some(cli::Commands::Map {
	    seq_files,
	    input_list,
	    index_prefix,
	    verbose,
        }) => {
	    init_log(if *verbose { 2 } else { 1 });
	    info!("Loading SBWT index...");

	    let (sbwt, lcs) = map::load_sbwt(index_prefix.clone().unwrap());

	    let translate_params = map::TranslateParams {
		k: match sbwt {
		    SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
			sbwt.n_kmers()
		    }
		},
		threshold: 14,
	    };

	    info!("Querying SBWT index...");
	    // TODO handle multiple files and `input_list`
	    let ms = map::query_sbwt(&seq_files[0], &sbwt, &lcs);

	    info!("Translating result...");
	    let ms_vec = ms.iter().map(|x| x.0).collect::<Vec<usize>>();
	    let runs = map::derandomize_ms(&ms_vec, &Some(translate_params.clone()));
	    let aln = map::translate_runs(&ms_vec, &runs, &Some(translate_params));
	    let run_lengths = map::run_lengths(&aln);

	    println!("q.start\tq.end\tlength\tmismatches");
	    run_lengths.iter().for_each(|x| println!("{}\t{}\t{}\t{}", x.0, x.1, x.2 + x.3, x.3));
	},
	None => {}
    }
}
