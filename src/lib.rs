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

//! kbo is an approximate local aligner based on converting [_k_-bounded matching
//! statistics](https://www.biorxiv.org/content/10.1101/2024.02.19.580943v1)
//! into a character representation of the underlying alignment sequence.
//!
//! Currently, kbo supports two main operations:
//!
//! - `kbo find` [matches](matches()) the _k_-mers in a query sequence with the
//! reference and reports the local alignment segments found within the
//! reference. Find is useful for problems that can be solved with
//! [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi).
//! - `kbo map` [maps](map()) the query sequence against a reference
//! sequence, and reports the nucleotide sequence of the alignment relative to
//! the reference. Map solves the same problem as
//! [snippy](https://github.com/tseemann/snippy) and [ska
//! map](https://docs.rs/ska/latest/ska/#ska-map).
//!
//! kbo uses the [Spectral Burrows-Wheeler
//! Transform](https://docs.rs/sbwt/latest/sbwt/) data structure that allows
//! efficient _k_-mer matching between a target and a query sequence and
//! fast retrieval of the _k_-bounded matching statistic for each _k_-mer match.
//!
//! # Installing the kbo executable
//! Run `cargo build --features cli --release`.
//!
//! # Usage
//!
//! kbo can be run directly on fasta files without an initial indexing step.
//! Prebuilt indexes are supported via `kbo build` but are only
//! relevant in `kbo find` analyses where the reference _k_-mers can be
//! concatenated into a single contig.
//!
//! ## kbo find
//!
//! To set up the example, download the fasta sequence of the [_Escherichia
//! coli_ Nissle
//! 1917](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000714595.1/) genome
//! from the NCBI and the [pks
//! island](https://raw.githubusercontent.com/tmaklin/clbtype/refs/heads/main/db/db.fasta)
//! gene sequences from GitHub. Example output was generated with versions
//! ASM71459v1 and rev 43bbd36.
//!
//! ### Find gene sequence locations
//! In the directory containing the input files, run
//! ```text
//! kbo find --reference db.fasta GCF_000714595.1_ASM71459v1_genomic.fna
//! ```
//! This will produce the output
//! 
//! ```text
//! TODO add output
//! ```
//!
//! ### Find presence/absence of gene sequences
//! Alternatively, if you are only interested in containment of the `db.fasta` genes in the assembly, run
//! ```text
//! kbo find --reference GCF_000714595.1_ASM71459v1_genomic.fna db.fasta
//! ```
//! which will return
//! ```text
//! TODO add output
//! ```text
//! ## kbo map
//! TODO write
//! ```
//!
//! ```
//!

#![warn(missing_docs,
        missing_debug_implementations, missing_copy_implementations,
        trivial_casts, trivial_numeric_casts,
        unsafe_code,
        unstable_features,
        unused_import_braces, unused_qualifications)]

use sbwt::SbwtIndexVariant;

pub mod derandomize;
pub mod format;
pub mod index;
pub mod translate;

/// Options and parameters for [Find](find)
#[non_exhaustive]
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FindOpts {
    /// Prefix match lengths with probability higher than `max_error_prob` to
    /// happen at random are considered noise.
    pub max_error_prob: f64,
    /// Maximum number of gap segments (gap opens) allowed before splitting an
    /// alignment.
    pub max_gaps: usize,
    /// Maximum length of a single gap segment before splitting an alignment.
    pub max_gap_len: usize,
}

impl Default for FindOpts {
    /// Default to these values:
    /// ```rust
    /// let mut opts = kbo::FindOpts::default();
    /// opts.max_error_prob = 0.0000001;
    /// opts.max_gaps = 0;
    /// opts.max_gap_len = 0;
    /// # let expected = kbo::FindOpts::default();
    /// # assert_eq!(opts, expected);
    /// ```
    ///
    fn default() -> FindOpts {
        FindOpts {
            max_error_prob: 0.0000001,
            max_gaps: 0,
            max_gap_len: 0,
        }
    }
}

/// Builds an SBWT index from some fasta or fastq files.
///
/// Reads all sequence data in `seq_files` and builds an SBWT index
/// with the parameters and resources specified in `build_opts` (see
/// [index::BuildOpts] for details).
///
/// Prebuilt indexes can currently only be used with kbo find.
///
/// All files and sequence data in `seq_files` are merged into the
/// same index. It is not possible extract the individual sequences
/// from the index after it has been built; run `kbo map -r
/// <query_file> <seq_files>` if you need to know which reference
/// sequences the alignments are for.
///
/// Returns a tuple containing the built
/// [SbwtIndexVariant](https://docs.rs/sbwt/latest/sbwt/enum.SbwtIndexVariant.html)
/// and
/// [sbwt::LcsArray](https://docs.rs/sbwt/latest/sbwt/struct.LcsArray.html).
///
/// Panics if a file in `seq_files` is not readable or a valid FASTX
/// file.
///
/// # Examples
/// ```rust
/// use kbo::build;
/// use kbo::index::BuildOpts;
///
/// let inputs: Vec<Vec<u8>> = vec![vec![b'A',b'A',b'A',b'G',b'A',b'A',b'C',b'C',b'A',b'-',b'T',b'C',b'A',b'G',b'G',b'G',b'C',b'G']];
///
/// let opts = BuildOpts::default();
/// let (sbwt_index, lcs_array) = build(&inputs, opts);
/// ```
///
pub fn build(
    seq_data: &[Vec<u8>],
    build_opts: index::BuildOpts,
) -> (SbwtIndexVariant, sbwt::LcsArray) {
    index::build_sbwt_from_vecs(seq_data, &Some(build_opts))
}

/// Matches a query fasta or fastq file against an SBWT index.
///
/// Queries the sequence data in `query_seq` against the SBWT index
/// `sbwt` and its LCS array `lcs` using [index::query_sbwt]. Then,
/// derandomizes the resulting _k_-bounded matching statistics vector
/// using [derandomize::derandomize_ms_vec] and translates the
/// matching statistics to a character representation of the alignment
/// using [translate::translate_ms_vec].
///
/// Returns a vector containing the character representation of the
/// alignment.
///
/// Panics if the query file is not readable or if it's not a valid
/// FASTX file.
///
/// # Output format
/// See the documentation for [translate].
///
/// # Example
/// ```rust
/// use kbo::build;
/// use kbo::matches;
/// use kbo::index::BuildOpts;
/// use kbo::derandomize::DerandomizeOpts;
///
/// let reference: Vec<Vec<u8>> = vec![vec![b'A',b'A',b'A',b'G',b'A',b'A',b'C',b'C',b'A',b'-',b'T',b'C',b'A',b'G',b'G',b'G',b'C',b'G']];
/// let mut opts = BuildOpts::default();
/// opts.k = 3;
/// let (sbwt, lcs) = build(&reference, opts);
///
/// let query = vec![b'G',b'T',b'G',b'A',b'C',b'T',b'A',b'T',b'G',b'A',b'G',b'G',b'A',b'T'];
///
/// let ms_vectors = matches(&query, &sbwt, &lcs, DerandomizeOpts::default());
/// // `ms_vectors` has ['-','-','-','-','-','-','-','-','-','M','M','M','-','-']
/// # assert_eq!(ms_vectors, vec!['-','-','-','-','-','-','-','-','-','M','M','M','-','-']);
/// ```
///
pub fn matches(
    query_seq: &[u8],
    sbwt: &SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
    derand_opts: derandomize::DerandomizeOpts,
) -> Vec<char> {
    let (k, threshold) = match sbwt {
	SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
	    (sbwt.k(), derandomize::random_match_threshold(sbwt.k(), sbwt.n_kmers(), 4_usize, derand_opts.max_error_prob))
	},
    };

    let noisy_ms: Vec<usize> = index::query_sbwt(query_seq, sbwt, lcs).iter().map(|x| x.0).collect();
    let derand_ms = derandomize::derandomize_ms_vec(&noisy_ms, k, threshold);

    translate::translate_ms_vec(&derand_ms, k, threshold)
}

/// Maps a query sequence against a reference sequence.
///
/// Maps the sequence data in `ref_seq` against the SBWT index
/// `query_sbwt` and `query_lcs` and converts the alignment to a
/// mapping relative to `ref_seq`.
///
/// Return the reference sequence with characters that are not present
/// in the query masked with a '-'.
///
/// # Examples
/// ```rust
/// use kbo::build;
/// use kbo::map;
/// use kbo::index::BuildOpts;
/// use kbo::derandomize::DerandomizeOpts;
///
/// let query: Vec<Vec<u8>> = vec![vec![b'A',b'A',b'A',b'G',b'A',b'A',b'C',b'C',b'A',b'-',b'T',b'C',b'A',b'G',b'G',b'G',b'C',b'G']];
/// let mut opts = BuildOpts::default();
/// opts.k = 3;
/// opts.build_select = true;
/// let (sbwt_query, lcs_query) = build(&query, opts);
///
/// let reference = vec![b'G',b'T',b'G',b'A',b'C',b'T',b'A',b'T',b'G',b'A',b'G',b'G',b'A',b'T'];
///
/// let alignment = map(&reference, &sbwt_query, &lcs_query, DerandomizeOpts::default());
/// // `ms_vectors` has [45,45,45,45,45,45,45,45,45,65,71,71,45,45]
/// # assert_eq!(alignment, vec![45,45,45,45,45,45,45,45,45,65,71,71,45,45]);
/// ```
///
pub fn map(
    ref_seq: &[u8],
    query_sbwt: &SbwtIndexVariant,
    query_lcs: &sbwt::LcsArray,
    derand_opts: derandomize::DerandomizeOpts,
) -> Vec<u8> {
    let (k, threshold) = match query_sbwt {
	SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
	    (sbwt.k(), derandomize::random_match_threshold(sbwt.k(), sbwt.n_kmers(), 4_usize, derand_opts.max_error_prob))
	},
    };

    let noisy_ms = index::query_sbwt(ref_seq, query_sbwt, query_lcs);
    let derand_ms = derandomize::derandomize_ms_vec(&noisy_ms.iter().map(|x| x.0).collect::<Vec<usize>>(), k, threshold);

    let translation = translate::translate_ms_vec(&derand_ms, k, threshold);
    let refined = translate::refine_translation(&translation, &noisy_ms, query_sbwt, threshold);

    format::relative_to_ref(ref_seq, &refined)
}

/// Finds the _k_-mers from an SBWT index in a query fasta or fastq file.
///
/// Aligns the sequence data in `query_seq` against the SBWT index
/// `sbwt` and its LCS array `lcs` using [matches][matches()]. Then uses
/// [format::run_lengths] to extract the local alignments from the
/// matching statistics.
///
/// Returns a vector of tuples, where each element represents a local
/// alignment block and contains the following values:
/// 1. Start of local alignment block in query (1-based indexing).
/// 2. End of local alignment block in query.
/// 3. Number of matches in the block.
/// 4. Number of mismatches and 1-character insertions in the block.
///
/// # Examples
///
/// TODO Add better examples to find()
///
/// ```rust
/// use kbo::build;
/// use kbo::find;
/// use kbo::FindOpts;
/// use kbo::index::BuildOpts;
/// use kbo::format::RLE;
///
/// let reference: Vec<Vec<u8>> = vec![vec![b'A',b'A',b'A',b'G',b'A',b'A',b'C',b'C',b'A',b'-',b'T',b'C',b'A',b'G',b'G',b'G',b'C',b'G']];
/// let mut opts = BuildOpts::default();
/// opts.k = 3;
/// let (sbwt, lcs) = build(&reference, opts);
///
/// let query = vec![b'G',b'T',b'G',b'A',b'C',b'T',b'A',b'T',b'G',b'A',b'G',b'G',b'A',b'T'];
///
/// let local_alignments = find(&query, &sbwt, &lcs, FindOpts::default());
/// // `local_alignments` has [(10, 12, 3, 0)]
/// # assert_eq!(local_alignments, vec![RLE{start: 10, end: 12, matches: 3, mismatches: 0, jumps: 0, gap_bases: 0, gap_opens: 0}]);
/// ```
///
pub fn find(
    query_seq: &[u8],
    sbwt: &SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
    find_opts: FindOpts,
) -> Vec<format::RLE> {
    let mut derand_opts = derandomize::DerandomizeOpts::default();
    derand_opts.max_error_prob = find_opts.max_error_prob.clone();
    let aln = matches(query_seq, sbwt, lcs, derand_opts);
    if find_opts.max_gaps > 0 || find_opts.max_gap_len > 0 {
        format::run_lengths_gapped(&aln, find_opts.max_gaps, find_opts.max_gap_len)
    } else {
        format::run_lengths(&aln)
    }
}
