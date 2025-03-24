# kbo
Local alignment search with _k_-bounded matching statistics.

kbo is an approximate local aligner based on converting [_k_-bounded matching
statistics](https://www.biorxiv.org/content/10.1101/2024.02.19.580943v1)
into a character representation of the underlying alignment sequence.

Documentation is available at [https://docs.rs/kbo](https://docs.rs/kbo).

## Usage
kbo is distributed as three separate Rust packages:
- [kbo](https://github.com/tmaklin/kbo) contains a Rust library implementing the core algorithm (this repository).
- [kbo-cli](https://github.com/tmaklin/kbo-cli) provides a command-line interface for `kbo call`, `kbo find` and `kbo map`.
- [kbo-gui](https://github.com/tmaklin/kbo-gui) is a work in progress WebAssembly graphical user interface for running kbo in the browser.

## About
kbo supports three main operations:

- `kbo call` calls single and multi base substitutions,
  insertions, and deletions in a query sequence against a reference and
  reports their positions and sequences. Call is useful for problems that
  require [.vcf files](https://samtools.github.io/hts-specs/VCFv4.2.pdf).
- `kbo find` matches the _k_-mers in a query sequence with the
  reference and reports the local alignment segments found within the
  reference. Find is useful for problems that can be solved with
  [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi).
- `kbo map` maps the query sequence against a reference
  sequence, and reports the nucleotide sequence of the alignment relative to
  the reference. Map solves the same problem as
  [snippy](https://github.com/tseemann/snippy) and [ska
  map](https://docs.rs/ska/latest/ska/#ska-map).

In addition to the three main operations, the core library provides high and
low-level functions that can be used in other libraries.

kbo uses the [Spectral Burrows-Wheeler
Transform](https://docs.rs/sbwt/latest/sbwt/) data structure that allows
efficient _k_-mer matching between a target and a query sequence and
fast retrieval of the _k_-bounded matching statistic for each _k_-mer match.

## License
kbo is dual-licensed under the [MIT](LICENSE-MIT) and [Apache 2.0](LICENSE-APACHE) licenses.
