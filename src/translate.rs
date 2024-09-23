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
//! Translating deterministic _k_-bounded matching statistics into alignments.
//!
//! The translated alignment is encoded using the following characters:
//! - **M** : Match between query and reference.
//! - **-** : Characters in the query that are not found in the reference.
//! - **X** : Single character mismatch or insertion into the query.
//! - **R** : Two consecutive 'R's signify a discontinuity in the alignment.
//!           The right 'R' is at the start of a _k_-mer that is not adjacent
//!           to the last character in the _k_-mer corresponding to the left
//!           'R'. This implies either a deletion of unknown length in the query,
//!           or insertion of _k_-mers from elsewhere in the reference into the query.
//!
//! ## Translation algorithm for _k_-bounded matching statistics
//! TODO Describe how the different MS vectors translate into alignments.

/// Translates a single derandomized _k_-bounded matching statistic.
///
/// Translates the current derandomized matching statistic (MS)
/// `ms_curr` based on the values of its left `ms_prev` and right
/// `ms_next` neighbors and the lower bound `threshold` for random
/// matches.
///
/// Returns a tuple containing the translation of the current MS and
/// translation of the right neighbor match. The right neighbor is an
/// empty string literal ' ' if translation of the current MS does not
/// affect its value.
///
/// # Examples
/// TODO add examples to translate_ms_val
pub fn translate_ms_val(
    ms_curr: i64,
    ms_next: i64,
    ms_prev: i64,
    threshold: usize,
) -> (char, char) {
    assert!(threshold > 1);

    let aln_curr: char;
    let mut aln_next: char = ' ';
    if ms_curr > threshold as i64 && ms_next > 0 && ms_next < threshold as i64 {
	// Current position is first character in a jump to another k-mer,
	// or there is deletion of unknown length in the query wrt. the reference.
	//
	// Use two consecutive 'R's to denote breakpoint between two k-mers
	aln_curr = 'R';
	aln_next = 'R';
    } else if ms_curr <= 0 {
	// Start of a mismatch region
	if ms_next == 1 && ms_prev > 0 {
	    // Mismatched character or insertion of 1 character in the query.
	    //
	    // Use 'X' for mismatch or 1 character insert
	    aln_curr = 'X';
	} else {
	    // Insertion of more than 1 characters in the query
	    //
	    // Use '-' to denote inserts of more than 1 characters
	    aln_curr = '-';
	}
    } else {
	// Other values are always a match, use 'M' for these
	aln_curr = 'M';
    }

    (aln_curr, aln_next)
}

/// Translates a sequence of derandomized _k_-bounded matching statistics.
///
/// Iterates over a derandomized sequence of _k_bounded matching
/// statistics `derand_ms` for _k_-mers with size `k` derandomized
/// with the threshold `threshold`.
///
/// Returns a sequence containing a character representation of the
/// underlying alignment.
///
/// # Examples
/// TODO Add examples to translate_ms_vec documentation.
///
pub fn translate_ms_vec(
    derand_ms: &[i64],
    k: usize,
    threshold: usize,
) -> Vec<char> {
    assert!(k > 0);
    assert!(threshold > 1);
    assert!(derand_ms.len() > 2);

    let len = derand_ms.len();
    let mut res = vec![' '; len];

    // Traverse the derandomized matchibng statistics
    for pos in 0..len {
	let prev: i64 = if pos > 1 { derand_ms[pos - 1] } else { 31 };
	let curr: i64 = derand_ms[pos];
	let next: i64 = if pos < len - 1 { derand_ms[pos + 1] } else { derand_ms[pos] };

	let (aln_curr, aln_next) = translate_ms_val(curr, next, prev, threshold);

	res[pos] = aln_curr;
	if pos + 1 < len - 1 && aln_next != ' ' {
	    res[pos + 1] = aln_next;
	}
    }

    res
}

////////////////////////////////////////////////////////////////////////////////
// Tests
//
#[cfg(test)]
mod tests {
    // Test cases for translate_ms_val
    // Comments use '-' for characters that are not in a ref or query sequence
    #[test]
    fn translate_ms_val_with_discontinuity() {
	// Parameters       : k = 3, threshold = 2
	//
	// Ref sequence     : A,C,G,T,T,T,C,A,G
	// Query sequence   : A,C,G,-,-,-,C,A,G
	//
	// Result MS vector : 1,2,3,1,2,3
	// Testing this pos :     |
	// Expected output  : M,M,R,R,M,M

	let expected = ('R','R');
	let got = super::translate_ms_val(3, 1, 2, 2);

	assert_eq!(got, expected);
    }

    #[test]
    fn translate_ms_val_with_mismatch() {
	// Parameters       : k = 3, threshold = 2
	//
	// Ref sequence     : A,C,G,T,C,A,G
	// Query sequence   : A,C,G,C,C,A,G
	//
	// Result MS vector : 1,2,3,0,1,2,3
	// Testing this pos :       |
	// Expected output  : M,M,M,X,M,M,M

	let expected = ('X',' ');
	let got = super::translate_ms_val(0, 1, 3, 2);

	assert_eq!(got, expected);
    }

    #[test]
    fn translate_ms_val_with_single_insertion() {
	// Parameters       : k = 3, threshold = 2
	//
	// Ref sequence     : A,C,G,-,C,A,G
	// Query sequence   : A,C,G,C,C,A,G
	//
	// Result MS vector : 1,2,3,0,1,2,3
	// Testing this pos :       |
	// Expected output  : M,M,M,X,M,M,M

	// Note this is identical to translate_ms_with_mismatch, these
	// are indistinguishible in outputs but have different inputs.
	// Kept here as a demonstration.
	let expected = ('X', ' ');
	let got = super::translate_ms_val(0, 1, 3, 2);

	assert_eq!(got, expected);
    }

    #[test]
    fn translate_ms_val_with_many_insertions() {
	// Parameters       : k = 3, threshold = 2
	//
	// Ref sequence     : A,C,G, -,-,C,A,G
	// Query sequence   : A,C,G, T,T,C,C,A,G
	//
	// Result MS vector : 1,2,3,-1,0,1,2,3
	// Testing this pos :        |
	// Expected output  : M,M,M, -,-,M,M,M

	let expected = ('-', ' ');
	let got = super::translate_ms_val(-1, 0, 3, 2);

	assert_eq!(got, expected);
    }

    #[test]
    fn translate_ms_val_with_only_matches() {
	// Parameters       : k = 3, threshold = 2
	//
	// Ref sequence     : A,C,G,C,A,G
	// Query sequence   : A,C,G,C,A,G
	//
	// Result MS vector : 1,2,3,3,3,3
	// Testing this pos :     |
	// Expected output  : M,M,M,M,M,M

	let expected = ('M', ' ');
	let got = super::translate_ms_val(1, 2, 3, 2);

	assert_eq!(got, expected);
    }
    
    #[test]
    fn translate_ms_vec() {
	// Parameters       : k = 3, threshold = 2
	// TODO check the k-mers
	// Ref sequence     : A,A,A,G,A,A,C,C,A,-,T,C,A, -,-,G,G,G,C,G
	// Query sequence   : C,A,A,G,-,-,C,C,A,C,T,C,A, T,T,G,G,G,T,C
	//
	// Result MS vector : 0,1,2,3,    1,2,3,0,1,2,3,-1,0,1,2,3,0,1
	// Expected output  : X,M,M,R,    R,M,M,X,M,M,M, -,-,M,M,M,-,-

	let input: Vec<i64> = vec![0,1,2,3,1,2,3,0,1,2,3,-1,0,1,2,3,0,1];
	let expected: Vec<char> = vec!['X','M','M','R','R','M','M','X','M','M','M','-','-','M','M','M','-','-'];
	let got = super::translate_ms_vec(&input, 3, 2);

	assert_eq!(got, expected);
    }
}
