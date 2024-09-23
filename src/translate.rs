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

    let mut aln_curr = ' ';
    let mut aln_next = ' ';
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
    for mut pos in 0..len {
	let prev: i64 = if pos > 1 { derand_ms[pos - 1] } else { 31 };
	let curr: i64 = derand_ms[pos];
	let next: i64 = if pos < len - 1 { derand_ms[pos + 1] } else { derand_ms[pos] };

	let mut aln_curr = res[pos];
	let mut aln_next = if pos + 1 < len - 1 { res[pos + 1] } else { 'M' };

	let (aln_curr, aln_next) = translate_ms_val(curr, next, prev, threshold);

	res[pos] = aln_curr;
	if pos + 1 < len - 1 && aln_next != ' ' {
	    res[pos + 1] = aln_next;
	}
    }

    return res;
}

////////////////////////////////////////////////////////////////////////////////
// Tests
//
#[cfg(test)]
mod tests {
    #[test]
    fn translate_ms_vec() {
	let expected = vec!['M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M'];
	// let input_ms: Vec<usize> = vec![1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,11,11,12,11,10,11,11,12,11,12,10,11,12,12,10,11,11,11,11,11,11,10,11,11,12,13,11,12,13,14,15,16,13,14,15,16,12,12,13,14,15,16,17,18,19,20,21,22,12,10,10,11,12,11,10,11,12,11,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,13,14,15,12,12,10,11,11,11,12,13,13,14,15,11,11,11,11,11,11,12,13,14,11,11,11,11,12,13,12,12,12,12,13,12,13,14,12,13,11,12,12,11,12,11,12,13,14,14,13,14,15,15,16,17,18,19,19,19,20,21,22,12,13,11,11,12,12,13,14,15,16,17,18,19,20,21,22,10,11,9,10,10,11,11,12,11,11,12,13,13,14,12,11,11,12,13,12,13,12,12,12,12,13,11,12,12,10,11,11,10,11,11,12,10,9,10,10,10,11,12,10,9,10,10,10,11,10,11,12,10,8,9,10,9,9,10,9,10,10,10,11,12,13,14,15,16,17,13,11,11,11,12,11,11,12,12,11,11,12,12,13,14,15,11,12,10,11,9,10,11,11,11,11,11,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,11,12,12,13,11,11,12,13,14,13,11,11,12,13,14,15,16,17,18,19,20,21,11,12,11,11,12,11,12,12,12,12,11,10,11,12,11,11,12,13,12,12,11,12,13,13,13,11,11,12,11,12,13,12,13,14,15,16,17,18,19,20,21,11,12,13,9,10,11,10,10,10,11,12,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27];
	let input_runs: Vec<i64> = vec![1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,-47,-46,-45,-44,-43,-42,-41,-40,-39,-38,-37,-36,-35,-34,-33,-32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,-164,-163,-162,-161,-160,-159,-158,-157,-156,-155,-154,-153,-152,-151,-150,-149,-148,-147,-146,-145,-144,-143,-142,-141,-140,-139,-138,-137,-136,-135,-134,-133,-132,-131,-130,-129,-128,-127,-126,-125,-124,-123,-122,-121,-120,-119,-118,-117,-116,-115,-114,-113,-112,-111,-110,-109,-108,-107,-106,-105,-104,-103,-102,-101,-100,-99,-98,-97,-96,-95,-94,-93,-92,-91,-90,-89,-88,-87,-86,-85,-84,-83,-82,-81,-80,-79,-78,-77,-76,-75,-74,-73,-72,-71,-70,-69,-68,-67,-66,-65,-64,-63,-62,-61,-60,-59,-58,-57,-56,-55,-54,-53,-52,-51,-50,-49,-48,-47,-46,-45,-44,-43,-42,-41,-40,-39,-38,-37,-36,-35,-34,-33,-32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,-63,-62,-61,-60,-59,-58,-57,-56,-55,-54,-53,-52,-51,-50,-49,-48,-47,-46,-45,-44,-43,-42,-41,-40,-39,-38,-37,-36,-35,-34,-33,-32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27];
	let got = super::translate_ms_vec(&input_runs, 31, 22);
	assert_eq!(got, expected);
    }
}
