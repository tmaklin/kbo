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
use std::io::Write;

use clap::Parser;
use log::info;
use needletail::Sequence;
use rayon::iter::ParallelIterator;
use rayon::iter::IntoParallelRefIterator;

// Command-line interface
mod cli;

/// Initializes the logger with verbosity given in `log_max_level`.
fn init_log(log_max_level: usize) {
    stderrlog::new()
	.module(module_path!())
	.quiet(false)
	.verbosity(log_max_level)
	.timestamp(stderrlog::Timestamp::Off)
	.init()
	.unwrap();
}

/// Use `sablast` to list the available commands or `sablast <command>` to run.
fn main() {
    let cli = cli::Cli::parse();

    // Subcommands:
    match &cli.command {
        Some(cli::Commands::Build {
	    seq_files,
	    output_prefix,
	    num_threads,
	    mem_gb,
	    temp_dir,
	    verbose,
        }) => {
	    init_log(if *verbose { 2 } else { 1 });

            let sbwt_build_options = sablast::index::BuildOpts {
		num_threads: *num_threads,
		mem_gb: *mem_gb,
		temp_dir: temp_dir.clone(),
                ..Default::default()
            };

	    info!("Building SBWT index from {} files...", seq_files.len());
	    let (sbwt, lcs) = sablast::build(seq_files, sbwt_build_options);

	    info!("Serializing SBWT index to {}.sbwt ...", output_prefix.as_ref().unwrap());
	    info!("Serializing LCS array to {}.lcs ...", output_prefix.as_ref().unwrap());
	    sablast::index::serialize_sbwt(output_prefix.as_ref().unwrap(), &sbwt, &lcs);

	},
        Some(cli::Commands::Find {
	    seq_files,
	    index_prefix,
	    num_threads,
	    verbose,
        }) => {
	    init_log(if *verbose { 2 } else { 1 });
	    rayon::ThreadPoolBuilder::new()
		.num_threads(*num_threads)
		.thread_name(|i| format!("rayon-thread-{}", i))
		.build_global()
		.unwrap();

	    info!("Loading SBWT index...");
	    let (sbwt, lcs) = sablast::index::load_sbwt(index_prefix.as_ref().unwrap());

	    info!("Querying SBWT index...");
	    println!("query\tref\tq.start\tq.end\tstrand\tlength\tmismatches");
	    seq_files.par_iter().for_each(|file| {

		let mut reader = needletail::parse_fastx_file(file).expect("valid path/file");
		let Some(rec) = reader.next() else { panic!("Invalid query {}", file); };
		let seqrec = rec.expect("Valid fastX record");
		let seq = seqrec.normalize(true);

		// Get local alignments
		let run_lengths = sablast::find(&seq, &sbwt, &lcs);

		// Print results with query and ref name added
		run_lengths.iter().for_each(|x| {
		    let stdout = std::io::stdout();
		    let _ = writeln!(&mut stdout.lock(),
				     "{}\t{}\t{}\t{}\t{}\t{}\t{}",
				     file, index_prefix.as_ref().unwrap(), x.0, x.1, x.2, x.3, x.4);
		});
	    });
	},
	None => {}
    }
}
