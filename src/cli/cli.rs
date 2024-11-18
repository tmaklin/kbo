// kbo: Spectral Burrows-Wheeler transform accelerated local alignment search
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
#[cfg(feature = "cli")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Subcommand)]
#[cfg(feature = "cli")]
pub enum Commands {
    // Build SBWT index
    Build {
        // Input fasta or fastq sequence file(s)
        #[arg(group = "input", required = true)]
        seq_files: Vec<String>,

        // Outputs
        #[arg(short = 'o', long = "output-prefix", required = false, help_heading = "Output")]
        output_prefix: Option<String>,

        // Build parameters
        // // k-mer size
        #[arg(short = 'k', default_value_t = 31, help_heading = "Build options")]
        kmer_size: usize,
        // // prefix precalc
        #[arg(short = 'p', long = "prefix-precalc", default_value_t = 8, help_heading = "Build options")]
        prefix_precalc: usize,
        // // deduplicate k-mer batches
        #[arg(short = 'd', long = "dedup-batches", default_value_t = false, help_heading = "Build options")]
        dedup_batches: bool,

        // Resources
        // // Threads
        #[arg(short = 't', long = "threads", default_value_t = 1)]
        num_threads: usize,
        // // Memory in GB
        #[arg(short = 'm', long = "mem-gb", default_value_t = 4, help_heading = "Build options")]
        mem_gb: usize,
        // // Temporary directory
        #[arg(long = "temp-dir", required = false, help_heading = "Build options")]
        temp_dir: Option<String>,

        // Verbosity
        #[arg(long = "verbose", default_value_t = false)]
        verbose: bool,
    },

    // Find indexed k-mers in a query
    Find {
        // Input fasta or fastq query file(s)
        #[arg(group = "input", required = true)]
        query_files: Vec<String>,

        // Reference
        // // Sequence file
        #[arg(short = 'r', long = "reference", group = "reference", help_heading = "Input")]
        ref_file: Option<String>,
        // // ... or a prebuilt index
        #[arg(short = 'i', long = "index", group = "reference", help_heading = "Input")]
        index_prefix: Option<String>,
        // // Concatenate contigs in reference
        #[arg(long = "detailed", help_heading = "Input", default_value_t = false)]
        detailed: bool,

        // Parameters
        // // Minimum length to report an alignment
        #[arg(long = "min-len", default_value_t = 100, help_heading = "Algorithm")]
        min_len: u64,
        #[arg(long = "max-gap-len", default_value_t = 0, help_heading = "Algorithm")]
        max_gap_len: u64,

        // // Upper bound for random match probability
        #[arg(long = "max-error-prob", default_value_t = 0.0000001, help_heading = "Algorithm")]
        max_error_prob: f64,

        // Resources
        // // Threads
        #[arg(short = 't', long = "threads", default_value_t = 1)]
        num_threads: usize,

        // Build parameters
        // // k-mer size
        #[arg(short = 'k', default_value_t = 31, help_heading = "Build options")]
        kmer_size: usize,
        // // prefix precalc
        #[arg(short = 'p', long = "prefix-precalc", default_value_t = 8, help_heading = "Build options")]
        prefix_precalc: usize,
        // // deduplicate k-mer batches
        #[arg(short = 'd', long = "dedup-batches", default_value_t = false, help_heading = "Build options")]
        dedup_batches: bool,
        // // Memory in GB
        #[arg(short = 'm', long = "mem-gb", default_value_t = 4, help_heading = "Build options")]
        mem_gb: usize,
        // // Temporary directory
        #[arg(long = "temp-dir", required = false, help_heading = "Build options")]
        temp_dir: Option<String>,


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

        // Parameters
        // // Upper bound for random match probability
        #[arg(long = "max-error-prob", default_value_t = 0.0000001, help_heading = "Algorithm")]
        max_error_prob: f64,

        // Resources
        // // Threads
        #[arg(short = 't', long = "threads", default_value_t = 1)]
        num_threads: usize,

        // Build parameters
        // // k-mer size
        #[arg(short = 'k', default_value_t = 31, help_heading = "Build options")]
        kmer_size: usize,
        // // prefix precalc
        #[arg(short = 'p', long = "prefix-precalc", default_value_t = 8, help_heading = "Build options")]
        prefix_precalc: usize,
        // // deduplicate k-mer batches
        #[arg(short = 'd', long = "dedup-batches", default_value_t = false, help_heading = "Build options")]
        dedup_batches: bool,
        // // Memory in GB
        #[arg(short = 'm', long = "mem-gb", default_value_t = 4, help_heading = "Build options")]
        mem_gb: usize,
        // // Temporary directory
        #[arg(long = "temp-dir", required = false, help_heading = "Build options")]
        temp_dir: Option<String>,

        // Verbosity
        #[arg(long = "verbose", default_value_t = false)]
        verbose: bool,
    },
}
