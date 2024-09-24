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
use clap::Parser;
use log::info;

// Command-line interface
mod cli;

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

            let sbwt_build_options = sablast::index::BuildOpts {
		num_threads: *num_threads,
		mem_gb: *mem_gb,
		temp_dir: temp_dir.clone(),
                ..Default::default()
            };
	    // TODO Handle --input-list in sablast build

	    let mut inputs = seq_files.clone();
	    info!("Reading input files...");
	    let mut seqs: Vec<Vec<u8>> = Vec::new();
	    inputs.iter().for_each(|file| {
		let mut reader = needletail::parse_fastx_file(file).unwrap_or_else(|_| panic!("Expected valid fastX file at {}", file));
		loop {
		    let rec = reader.next();
		    match rec {
			Some(Ok(seqrec)) => {
			    seqs.push(seqrec.normalize(true).as_ref().to_vec());
			},
			_ => break
		    }
		}
	    });

	    info!("Building SBWT index...");
	    let (sbwt, lcs) = sablast::index::build_sbwt_from_vecs(&seqs, &Some(sbwt_build_options));

	    info!("Serializing SBWT index...");
	    sablast::index::serialize_sbwt(&output_prefix.as_ref().unwrap(), &sbwt, &lcs);

	},
        Some(cli::Commands::Map {
	    seq_files,
	    input_list,
	    index_prefix,
	    verbose,
        }) => {
	    init_log(if *verbose { 2 } else { 1 });
	    info!("Loading SBWT index...");

	    let (sbwt, lcs) = sablast::index::load_sbwt(index_prefix.as_ref().unwrap());

	    // TODO Handle `--input-list in sablast map

	    // TODO Query multiple inputs in sablast map
	    info!("Querying SBWT index...");

	    let aln = sablast::map(&seq_files[0], &sbwt, &lcs);
	    let mut run_lengths: Vec<(usize, usize, usize, usize, char)> = sablast::format::run_lengths(&aln.0).iter().map(|x| (x.0, x.1, x.2, x.3, '+')).collect();
	    let mut run_lengths_rev: Vec<(usize, usize, usize, usize, char)> = sablast::format::run_lengths(&aln.1).iter().map(|x| (x.0, x.1, x.2, x.3, '-')).collect();
	    run_lengths.append(&mut run_lengths_rev);

	    run_lengths.sort_by_key(|x| x.0);

	    println!("query\tref\tq.start\tq.end\tstrand\tlength\tmismatches");
	    run_lengths.iter().for_each(|x| println!("{}\t{}\t{}\t{}\t{}\t{}\t{}", &seq_files[0], &index_prefix.clone().unwrap(), x.0, x.1, x.4, x.2 + x.3, x.3));
	},
	None => {}
    }
}
