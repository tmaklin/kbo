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

// Reads all sequence data from a fastX file
fn read_fastx_file(
    file: &str,
) -> Vec<Vec<u8>> {
    let mut seq_data: Vec<Vec<u8>> = Vec::new();
    let mut reader = needletail::parse_fastx_file(file).unwrap_or_else(|_| panic!("Expected valid fastX file at {}", file));
    loop {
	let rec = reader.next();
	match rec {
	    Some(Ok(seqrec)) => {
		seq_data.push(seqrec.normalize(true).as_ref().to_vec());
	    },
	    _ => break
	}
    }
    seq_data
}

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
///
/// # Input format detection
/// The sequence data is read using
/// [needletail::parser::parse_fastx_file](https://docs.rs/needletail/latest/needletail/parser/fn.parse_fastx_file.html).
///
/// Input file format (fasta or fastq) is detected automatically and
/// the files may be compressed in a
/// [DEFLATE-based](https://en.wikipedia.org/wiki/Deflate) format (.gz
/// files).
///
/// [Bzip2](https://sourceware.org/bzip2/) and
/// [liblzma](https://tukaani.org/xz/) compression (.bz2 and .xz
/// files) can be enabled using the needletail features field in
/// sablast Cargo.toml if compiling from source.
///
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
	    let mut seq_data: Vec<Vec<u8>> = Vec::new();
	    seq_files.iter().for_each(|file| {
		seq_data.append(&mut read_fastx_file(file));
	    });

	    let (sbwt, lcs) = sablast::build(&seq_data, sbwt_build_options);

	    info!("Serializing SBWT index to {}.sbwt ...", output_prefix.as_ref().unwrap());
	    info!("Serializing LCS array to {}.lcs ...", output_prefix.as_ref().unwrap());
	    sablast::index::serialize_sbwt(output_prefix.as_ref().unwrap(), &sbwt, &lcs);

	},
        Some(cli::Commands::Find {
	    query_files,
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
	    println!("query\tref\tq.start\tq.end\tstrand\tlength\tmismatches\tin.contig");
	    let stdout = std::io::stdout();
	    query_files.par_iter().for_each(|file| {

		let mut reader = needletail::parse_fastx_file(file).expect("valid path/file");
		while let Some(rec) = reader.next() {
		    let seqrec = rec.expect("Valid fastX record");
		    let contig = seqrec.id();
		    let seq = seqrec.normalize(true);

		    // Get local alignments for forward strand
		    let mut run_lengths: Vec<(usize, usize, char, usize, usize)> = sablast::find(&seq, &sbwt, &lcs).iter().map(|x| (x.0, x.1, '+', x.2 + x.3, x.3)).collect();

		    // Add local alignments for reverse _complement
		    run_lengths.append(&mut sablast::find(&seq.reverse_complement(), &sbwt, &lcs).iter().map(|x| (x.0, x.1, '-', x.2 + x.3, x.3)).collect());

		    // Sort by q.start
		    run_lengths.sort_by_key(|x| x.0);

		    // Print results with query and ref name added
		    run_lengths.iter().for_each(|x| {
			let _ = writeln!(&mut stdout.lock(),
					 "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
					 file, index_prefix.as_ref().unwrap(), x.0, x.1, x.2, x.3, x.4, std::str::from_utf8(contig).expect("UTF-8"));
		    });
		}
	    });
	},
        Some(cli::Commands::Map {
	    query_files,
	    ref_file,
	    num_threads,
	    verbose,
        }) => {
	    init_log(if *verbose { 2 } else { 1 });
	    rayon::ThreadPoolBuilder::new()
		.num_threads(*num_threads)
		.thread_name(|i| format!("rayon-thread-{}", i))
		.build_global()
		.unwrap();

	    let ref_data = read_fastx_file(ref_file);
	    let opts = sablast::index::BuildOpts { add_revcomp: true, build_select: true, ..Default::default() };

	    let stdout = std::io::stdout();
	    query_files.par_iter().for_each(|query_file| {
		let query_data = read_fastx_file(query_file);
		let (sbwt, lcs) = sablast::index::build_sbwt_from_vecs(&query_data, &Some(opts.clone()));

		let mut res: Vec<u8> = Vec::new();
		ref_data.iter().for_each(|ref_contig| {
		    res.append(&mut sablast::map(ref_contig, &sbwt, &lcs));
		});

		let _ = writeln!(&mut stdout.lock(),
				 ">{}\n{}", query_file, std::str::from_utf8(&res).expect("UTF-8"));
	    });
	},
	None => {}
    }
}
