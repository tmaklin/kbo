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
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(version)]
#[command(propagate_version = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Subcommand)]
pub enum Commands {
    // Build SBWT index
    Build {
        // Input fasta or fastq sequence file(s)
        #[arg(group = "input", required = true)]
        seq_files: Vec<String>,

	// Outputs
        #[arg(short = 'o', long = "out-prefix", required = false, help_heading = "Output")]
        output_prefix: Option<String>,

        // Resources
	// // Threads
        #[arg(short = 't', long = "threads", default_value_t = 1)]
        num_threads: usize,
	// // Memory in GB
        #[arg(short = 'm', long = "memory", default_value_t = 4)]
        mem_gb: usize,
	// // Temporary directory
        #[arg(long = "tmp-dir", required = false)]
        temp_dir: Option<String>,

	// Verbosity
        #[arg(long = "verbose", default_value_t = false)]
        verbose: bool,
    },

    // Find indexed k-mers in a query
    Find {
	// Input fasta or fastq query file(s)
        #[arg(group = "input", required = true)]
        seq_files: Vec<String>,

	// Index name
        #[arg(short = 'i', long = "index", required = true, help_heading = "Index")]
        index_prefix: Option<String>,

	// Resources
	// // Threads
        #[arg(short = 't', long = "threads", default_value_t = 1)]
        num_threads: usize,

	// Verbosity
        #[arg(long = "verbose", default_value_t = false)]
        verbose: bool,
    },

    // Map a query or queries to a reference and return the alignment
    Map {
	// Input fasta or fastq query file(s)
        #[arg(group = "input", required = true)]
        query_files: Vec<String>,

	// Reference fasta
        #[arg(short = 'r', long = "reference", required = true, help_heading = "Input")]
        ref_file: String,

	// Resources
	// // Threads
        #[arg(short = 't', long = "threads", default_value_t = 1)]
        num_threads: usize,

	// Verbosity
        #[arg(long = "verbose", default_value_t = false)]
        verbose: bool,
    },
}
