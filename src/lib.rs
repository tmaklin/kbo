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
//! Currently, kbo supports three main operations:
//!
//! - `kbo call` [calls](call()) single and multi base substitutions,
//!   insertions, and deletions in a query sequence against a reference and
//!   reports their positions and sequences. Call is useful for problems that
//!   require [.vcf files](https://samtools.github.io/hts-specs/VCFv4.2.pdf).
//! - `kbo find` [matches](matches()) the _k_-mers in a query sequence with the
//!   reference and reports the local alignment segments found within the
//!   reference. Find is useful for problems that can be solved with
//!   [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi).
//! - `kbo map` [maps](map()) the query sequence against a reference
//!   sequence, and reports the nucleotide sequence of the alignment relative to
//!   the reference. Map solves the same problem as
//!   [snippy](https://github.com/tseemann/snippy) and [ska
//!   map](https://docs.rs/ska/latest/ska/#ska-map).
//!
//! kbo uses the [Spectral Burrows-Wheeler
//! Transform](https://docs.rs/sbwt/latest/sbwt/) data structure that allows
//! efficient _k_-mer matching between a target and a query sequence and
//! fast retrieval of the _k_-bounded matching statistic for each _k_-mer match.
//!
//! # Installing the kbo executable
//! See installation instructions at [GitHub](https://github.com/tmaklin/kbo-cli).
//!
//! # Usage
//!
//! kbo can be run directly on fasta files without an initial indexing step.
//! Prebuilt indexes are supported via `kbo build` but are only
//! relevant in `kbo find` analyses where the reference _k_-mers can be
//! concatenated into a single contig.
//!
//! kbo can read inputs compressed in the DEFLATE format (gzip, zlib, etc.).
//! bzip2 and xz support can be enabled by adding the "bzip2" and "xz" feature
//! flags to [needletail](https://docs.rs/needletail) in the kbo Cargo.toml.
//!
//! ## kbo call
//!
//! Set up the example by downloading the fasta file for the [_Streptococcus
//! pneumoniae_
//! Spn23F](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000026665.1/)
//! genome from the NCBI and the [_S. pneumoniae_
//! 6952_7#3](https://www.ebi.ac.uk/ena/browser/view/GCA_001156685.2) assembly
//! from the ENA.
//!
//! ### Calling variants in a reference genome
//! In the directory with the downloaded files, run
//! ```text
//! kbo call --reference GCF_000026665.1_ASM2666v1_genomic.fna GCA_001156685.2.fasta.gz > variants.vcf
//! ```
//! This will write the variants in the [vcf v4.4](https://samtools.github.io/hts-specs/VCFv4.4.pdf) format
//!
//! <details>
//! <summary>
//! (click to view the first 20 lines)
//! </summary>
//!
//! ```text
//! ##fileformat=VCFv4.4
//! ##contig=<ID=NC_011900.1,length=2221315>
//! ##fileDate=20250324
//! ##source=kbo-cli v0.1.1
//! ##reference=GCF_000026665.1_ASM2666v1_genomic.fna
//! ##phasing=none
//! #CHROM          POS     ID  REF  ALT  QUAL  FILTER  INFO   FORMAT  unknown
//! NC_011900.1     83      .   G    A    .     .       .      GT      1
//! NC_011900.1     845     .   A    C    .     .       .      GT      1
//! NC_011900.1     1064    .   G    A    .     .       .      GT      1
//! NC_011900.1     1981    .   G    A    .     .       .      GT      1
//! NC_011900.1     2392    .   C    T    .     .       .      GT      1
//! NC_011900.1     2746    .   C    T    .     .       .      GT      1
//! NC_011900.1     3236    .   T    C    .     .       .      GT      1
//! NC_011900.1     3397    .   A    G    .     .       .      GT      1
//! NC_011900.1     3993    .   C    T    .     .       .      GT      1
//! NC_011900.1     4335    .   AA   A    .     .       INDEL  GT      1
//! NC_011900.1     4504    .   C    A    .     .       .      GT      1
//! NC_011900.1     4861    .   A    G    .     .       .      GT      1
//! NC_011900.1     5007    .   A    T    .     .       .      GT      1
//! ```
//!
//! </details>
//!
//! ## kbo find
//!
//! First download the fasta sequence of the [_Escherichia
//! coli_ Nissle
//! 1917](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000714595.1/) genome
//! from the NCBI and the [pks
//! island](https://raw.githubusercontent.com/tmaklin/clbtype/021e09f07ba43f79cc9c2f365be4058cebbcd7ce/db/IHE3034_pks_island_genes.fasta)
//! gene sequences from GitHub. Example output was generated with versions
//! ASM71459v1 and rev 021e09f.
//!
//! ### Find gene sequence locations
//! In the directory containing the input files, run
//! ```text
//! kbo find --max-gap-len 100 --reference IHE3034_pks_island_genes.fasta GCF_000714595.1_ASM71459v1_genomic.fna
//! ```
//! <details>
//! <summary>
//! This will produce the output (click to expand)
//! </summary>
//!
//! |query|ref|q.start|q.end|strand|length|mismatches|gap_bases|gap_opens|identity|coverage|query.contig|ref.contig|
//! |-----|---|-------|-----|------|------|----------|---------|---------|--------|--------|------------|----------|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2289596|2290543|+|948|0|0|0|100.00|1.90|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|IHE3034_pks_island_genes.fasta|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2239798|2289162|-|49365|7|367|12|99.24|98.06|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|IHE3034_pks_island_genes.fasta|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|5145962|5149449|+|3488|0|61|1|98.25|6.86|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|IHE3034_pks_island_genes.fasta|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|5354674|5356713|+|2040|1|0|0|99.95|4.08|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|IHE3034_pks_island_genes.fasta|
//!
//! </details>
//!
//! ### Find gene sequence locations with names
//! If you need to know which gene in db.fasta the matches are for, add the `--detailed` toggle:
//! ```text
//! kbo find --detailed --reference IHE3034_pks_island_genes.fasta GCF_000714595.1_ASM71459v1_genomic.fna
//! ```
//!
//! <details>
//! <summary>
//! This replaces the query.contig column with the name of the contig (click to expand)
//! </summary>
//!
//! |query|ref|q.start|q.end|strand|length|mismatches|gap_bases|gap_opens|identity|coverage|query.contig|ref.contig|
//! |-----|---|-------|-----|------|------|----------|---------|---------|--------|--------|------------|----------|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2289596|2289808|+|213|0|0|0|100.00|100.00|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbR\|locus_tag=ECOK1_RS11410\|product="colibactin biosynthesis LuxR family transcriptional regulator ClbR"\|protein_id=WP_000357141.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2289809|2290543|+|735|0|0|0|100.00|100.00|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbA\|locus_tag=ECOK1_RS11415\|product="colibactin biosynthesis phosphopantetheinyl transferase ClbA"\|protein_id=WP_001217110.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2279541|2289162|-|9622|1|0|0|99.99|100.01|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbB\|locus_tag=ECOK1_RS11405\|product="colibactin hybrid non-ribosomal peptide synthetase/type I polyketide synthase ClbB"\|protein_id=WP_001518711.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2276900|2279500|-|2601|0|0|0|100.00|100.00|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbC\|locus_tag=ECOK1_RS11400\|product="colibactin polyketide synthase ClbC"\|protein_id=WP_001297908.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2276021|2276890|-|870|0|0|0|100.00|100.00|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbD\|locus_tag=ECOK1_RS11395\|product="colibactin biosynthesis dehydrogenase ClbD"\|protein_id=WP_000982270.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2275743|2275991|-|249|0|0|0|100.00|100.00|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbE\|locus_tag=ECOK1_RS11390\|product="colibactin biosynthesis aminomalonyl-acyl carrier protein ClbE"\|protein_id=WP_001297917.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2274609|2275739|-|1131|0|0|0|100.00|100.00|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbF\|locus_tag=ECOK1_RS11385\|product="colibactin biosynthesis dehydrogenase ClbF"\|protein_id=WP_000337350.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2273344|2274612|-|1269|1|0|0|99.92|100.00|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbG\|locus_tag=ECOK1_RS11380\|product="colibactin biosynthesis acyltransferase ClbG"\|protein_id=WP_000159201.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2268500|2273296|-|4797|2|0|0|99.96|100.00|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbH\|locus_tag=ECOK1_RS11375\|product="colibactin non-ribosomal peptide synthetase ClbH"\|protein_id=WP_001304254.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2265418|2268450|-|3033|0|0|0|100.00|100.00|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbI\|locus_tag=ECOK1_RS11370\|product="colibactin polyketide synthase ClbI"\|protein_id=WP_000829570.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2258874|2265374|-|6501|0|0|0|100.00|100.00|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbJ\|locus_tag=ECOK1_RS11365\|product="colibactin non-ribosomal peptide synthetase ClbJ"\|protein_id=WP_001468003.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2259498|2260784|-|1287|2|0|0|99.84|19.91|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbK\|locus_tag=ECOK1_RS11360\|product="colibactin hybrid non-ribosomal peptide synthetase/type I polyketide synthase ClbK"\|protein_id=WP_000222467.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2252399|2258863|-|6465|2|0|0|99.97|100.00|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbK\|locus_tag=ECOK1_RS11360\|product="colibactin hybrid non-ribosomal peptide synthetase/type I polyketide synthase ClbK"\|protein_id=WP_000222467.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2253845|2255131|-|1287|1|0|0|99.92|19.80|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbJ\|locus_tag=ECOK1_RS11365\|product="colibactin non-ribosomal peptide synthetase ClbJ"\|protein_id=WP_001468003.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2250943|2252406|-|1464|0|0|0|100.00|100.00|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbL\|locus_tag=ECOK1_RS11355\|product="colibactin biosynthesis amidase ClbL"\|protein_id=WP_001297937.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2249442|2250881|-|1440|0|0|0|100.00|100.00|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbM\|locus_tag=ECOK1_RS11350\|product="precolibactin export MATE transporter ClbM"\|protein_id=WP_000217768.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2245077|2249445|-|4369|1|0|0|99.98|100.02|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbN\|locus_tag=ECOK1_RS11345\|product="colibactin non-ribosomal peptide synthetase ClbN"\|protein_id=WP_001327259.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2242587|2245046|-|2460|0|0|0|100.00|100.00|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbO\|locus_tag=ECOK1_RS11340\|product="colibactin polyketide synthase ClbO"\|protein_id=WP_001029878.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2241060|2242574|-|1515|0|0|0|100.00|100.00|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbP\|locus_tag=ECOK1_RS11335\|product="precolibactin peptidase ClbP"\|protein_id=WP_002430641.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2240345|2241067|-|723|0|0|0|100.00|100.00|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbQ\|locus_tag=ECOK1_RS11330\|product="colibactin biosynthesis thioesterase ClbQ"\|protein_id=WP_000065646.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|2239798|2240310|-|513|1|0|0|99.81|100.00|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbS\|locus_tag=ECOK1_RS11325\|product="colibactin self-protection protein ClbS"\|protein_id=WP_000290498.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|5145962|5147210|+|1249|0|0|0|100.00|85.31|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbL\|locus_tag=ECOK1_RS11355\|product="colibactin biosynthesis amidase ClbL"\|protein_id=WP_001297937.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|5147272|5148479|+|1208|0|0|0|100.00|83.89|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbM\|locus_tag=ECOK1_RS11350\|product="precolibactin export MATE transporter ClbM"\|protein_id=WP_000217768.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|5148478|5149449|+|972|0|0|0|100.00|22.25|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbN\|locus_tag=ECOK1_RS11345\|product="colibactin non-ribosomal peptide synthetase ClbN"\|protein_id=WP_001327259.1|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|IHE3034_pks_island_genes.fasta|5354674|5356713|+|2040|1|0|0|99.95|46.70|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|clbN\|locus_tag=ECOK1_RS11345\|product="colibactin non-ribosomal peptide synthetase ClbN"\|protein_id=WP_001327259.1|
//!
//! </details>
//!
//! Note that the current implementation `--detailed` slows down
//! the algorithm. Future versions of kbo may address this by incorporating
//! colors in the index structure.
//!
//! ### Find containment of gene sequences in assembly
//!
//! Alternatively, if you are only interested in whether the contigs in
//! `db.fasta` are present in the assembly, swap the reference and query above
//! run
//!
//! ```text
//! kbo find --reference GCF_000714595.1_ASM71459v1_genomic.fna IHE3034_pks_island_genes.fasta
//! ```
//!
//! <details>
//! <summary>
//! which will return (click to expand)
//! </summary>
//!
//! |query|ref|q.start|q.end|strand|length|mismatches|gap_bases|gap_opens|identity|coverage|query.contig|ref.contig|
//! |-----|---|-------|-----|------|------|----------|---------|---------|--------|--------|------------|----------|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|513|-|513|1|0|0|99.81|0.01|clbS\|locus_tag=ECOK1_RS11325\|product="colibactin self-protection protein ClbS"\|protein_id=WP_000290498.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|723|-|723|0|0|0|100.00|0.01|clbQ\|locus_tag=ECOK1_RS11330\|product="colibactin biosynthesis thioesterase ClbQ"\|protein_id=WP_000065646.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|1515|-|1515|0|0|0|100.00|0.03|clbP\|locus_tag=ECOK1_RS11335\|product="precolibactin peptidase ClbP"\|protein_id=WP_002430641.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|2460|-|2460|0|0|0|100.00|0.05|clbO\|locus_tag=ECOK1_RS11340\|product="colibactin polyketide synthase ClbO"\|protein_id=WP_001029878.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|4368|-|4368|1|0|0|99.98|0.08|clbN\|locus_tag=ECOK1_RS11345\|product="colibactin non-ribosomal peptide synthetase ClbN"\|protein_id=WP_001327259.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|1208|+|1208|0|0|0|100.00|0.02|clbM\|locus_tag=ECOK1_RS11350\|product="precolibactin export MATE transporter ClbM"\|protein_id=WP_000217768.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|1440|-|1440|0|0|0|100.00|0.03|clbM\|locus_tag=ECOK1_RS11350\|product="precolibactin export MATE transporter ClbM"\|protein_id=WP_000217768.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|1464|-|1464|0|0|0|100.00|0.03|clbL\|locus_tag=ECOK1_RS11355\|product="colibactin biosynthesis amidase ClbL"\|protein_id=WP_001297937.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|6465|-|6465|2|0|0|99.97|0.12|clbK\|locus_tag=ECOK1_RS11360\|product="colibactin hybrid non-ribosomal peptide synthetase/type I polyketide synthase ClbK"\|protein_id=WP_000222467.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|6501|-|6501|0|0|0|100.00|0.12|clbJ\|locus_tag=ECOK1_RS11365\|product="colibactin non-ribosomal peptide synthetase ClbJ"\|protein_id=WP_001468003.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|3033|-|3033|0|0|0|100.00|0.06|clbI\|locus_tag=ECOK1_RS11370\|product="colibactin polyketide synthase ClbI"\|protein_id=WP_000829570.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|4797|-|4797|2|0|0|99.96|0.09|clbH\|locus_tag=ECOK1_RS11375\|product="colibactin non-ribosomal peptide synthetase ClbH"\|protein_id=WP_001304254.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|1269|-|1269|1|0|0|99.92|0.02|clbG\|locus_tag=ECOK1_RS11380\|product="colibactin biosynthesis acyltransferase ClbG"\|protein_id=WP_000159201.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|1131|-|1131|0|0|0|100.00|0.02|clbF\|locus_tag=ECOK1_RS11385\|product="colibactin biosynthesis dehydrogenase ClbF"\|protein_id=WP_000337350.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|249|-|249|0|0|0|100.00|0.00|clbE\|locus_tag=ECOK1_RS11390\|product="colibactin biosynthesis aminomalonyl-acyl carrier protein ClbE"\|protein_id=WP_001297917.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|870|-|870|0|0|0|100.00|0.02|clbD\|locus_tag=ECOK1_RS11395\|product="colibactin biosynthesis dehydrogenase ClbD"\|protein_id=WP_000982270.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|2601|-|2601|0|0|0|100.00|0.05|clbC\|locus_tag=ECOK1_RS11400\|product="colibactin polyketide synthase ClbC"\|protein_id=WP_001297908.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|9621|-|9621|1|0|0|99.99|0.18|clbB\|locus_tag=ECOK1_RS11405\|product="colibactin hybrid non-ribosomal peptide synthetase/type I polyketide synthase ClbB"\|protein_id=WP_001518711.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|213|+|213|0|0|0|100.00|0.00|clbR\|locus_tag=ECOK1_RS11410\|product="colibactin biosynthesis LuxR family transcriptional regulator ClbR"\|protein_id=WP_000357141.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|735|+|735|0|0|0|100.00|0.01|clbA\|locus_tag=ECOK1_RS11415\|product="colibactin biosynthesis phosphopantetheinyl transferase ClbA"\|protein_id=WP_001217110.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|216|1464|+|1249|0|0|0|100.00|0.02|clbL\|locus_tag=ECOK1_RS11355\|product="colibactin biosynthesis amidase ClbL"\|protein_id=WP_001297937.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//! |IHE3034_pks_island_genes.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1156|4167|+|3012|1|0|0|99.97|0.06|clbN\|locus_tag=ECOK1_RS11345\|product="colibactin non-ribosomal peptide synthetase ClbN"\|protein_id=WP_001327259.1|GCF_000714595.1_ASM71459v1_genomic.fna|
//!
//! </details>
//!
//! ## kbo map
//!
//! kbo map can be used to align a query sequence against a reference sequence.
//! This is useful in for example generating a reference-based alignment of
//! multiple related genomes against a good reference assembly.
//!
//! To run this example, download the genome sequence of the [_E. coli_
//! UTI89](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000013265.1/) strain
//! from the NCBI (ASM1326v1) and [_E. coli_ Nissle
//! 1917](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000714595.1/)
//! (ASM71459v1).
//!
//! ### Reference-based alignment
//! Run
//! ```text
//! kbo map --reference GCF_000714595.1_ASM71459v1_genomic.fna GCF_000013265.1_ASM1326v1_genomic.fna > result.aln
//! ```
//!
//! which will write the alignment sequence to `result.aln`.
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
pub mod gap_filling;
pub mod index;
pub mod translate;
pub mod variant_calling;

/// Options and parameters for [call]
#[non_exhaustive]
#[derive(Clone, Debug)]
pub struct CallOpts {
    /// Prefix match lengths with probability higher than `max_error_prob` to
    /// happen at random are considered noise.
    pub max_error_prob: f64,

    /// [Build options](index::BuildOpts) for indexing the reference sequence, `k` must
    /// match the _k_-mer size used in indexing the query.
    pub sbwt_build_opts: index::BuildOpts,
}

impl Default for CallOpts {
    /// Default to these values:
    /// ```rust
    /// let mut opts = kbo::CallOpts::default();
    /// opts.max_error_prob = 0.0000001;
    /// opts.sbwt_build_opts = kbo::index::BuildOpts::default();
    /// opts.sbwt_build_opts.build_select = true;
    /// # let expected = kbo::CallOpts::default();
    /// # assert_eq!(opts.max_error_prob, expected.max_error_prob);
    /// # assert_eq!(opts.sbwt_build_opts.k, expected.sbwt_build_opts.k);
    /// # assert_eq!(opts.sbwt_build_opts.add_revcomp, expected.sbwt_build_opts.add_revcomp);
    /// # assert_eq!(opts.sbwt_build_opts.num_threads, expected.sbwt_build_opts.num_threads);
    /// # assert_eq!(opts.sbwt_build_opts.prefix_precalc, expected.sbwt_build_opts.prefix_precalc);
    /// # assert_eq!(opts.sbwt_build_opts.build_select, true);
    /// # assert_eq!(opts.sbwt_build_opts.mem_gb, expected.sbwt_build_opts.mem_gb);
    /// # assert_eq!(opts.sbwt_build_opts.dedup_batches, expected.sbwt_build_opts.dedup_batches);
    /// # assert_eq!(opts.sbwt_build_opts.temp_dir, expected.sbwt_build_opts.temp_dir);
    /// ```
    ///
    fn default() -> CallOpts {
        CallOpts {
            max_error_prob: 0.0000001,
            sbwt_build_opts: index::BuildOpts { build_select: true, ..Default::default() },
        }
    }
}

/// Options and parameters for [find]
#[non_exhaustive]
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FindOpts {
    /// Prefix match lengths with probability higher than `max_error_prob` to
    /// happen at random are considered noise.
    pub max_error_prob: f64,
    /// Maximum length of a single gap segment before splitting an alignment.
    pub max_gap_len: usize,
}

impl Default for FindOpts {
    /// Default to these values:
    /// ```rust
    /// let mut opts = kbo::FindOpts::default();
    /// opts.max_error_prob = 0.0000001;
    /// opts.max_gap_len = 0;
    /// # let expected = kbo::FindOpts::default();
    /// # assert_eq!(opts, expected);
    /// ```
    ///
    fn default() -> FindOpts {
        FindOpts {
            max_error_prob: 0.0000001,
            max_gap_len: 0,
        }
    }
}

/// Options and parameters for [matches](matches())
#[non_exhaustive]
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct MatchOpts {
    /// Prefix match lengths with probability higher than `max_error_prob` to
    /// happen at random are considered noise.
    pub max_error_prob: f64,
}

impl Default for MatchOpts {
    /// Default to these values:
    /// ```rust
    /// let mut opts = kbo::MatchOpts::default();
    /// opts.max_error_prob = 0.0000001;
    /// # let expected = kbo::MatchOpts::default();
    /// # assert_eq!(opts, expected);
    /// ```
    ///
    fn default() -> MatchOpts {
        MatchOpts {
            max_error_prob: 0.0000001,
        }
    }
}

/// Options and parameters for [map]
#[non_exhaustive]
#[derive(Clone, Debug)]
pub struct MapOpts {
    /// Prefix match lengths with probability higher than `max_error_prob` to
    /// happen at random are considered noise.
    pub max_error_prob: f64,

    /// Attempt to resolve gaps in the alignment with [gap
    /// filling](gap_filling::fill_gaps).
    pub fill_gaps: bool,

    /// Resolve substitutions, insertions, and deletions with [variant
    /// calling](call).
    pub call_variants: bool,

    /// Replace characters used by the internal representation of [translate]
    /// with nucleotide codes and gaps.
    pub format: bool,

    /// [Build options](index::BuildOpts) for SBWT, used if variant calling is
    /// requested. `k` must match the _k_-mer size used in indexing the query.
    pub sbwt_build_opts: index::BuildOpts,
}

impl Default for MapOpts {
    /// Default to these values:
    /// ```rust
    /// let mut opts = kbo::MapOpts::default();
    /// opts.max_error_prob = 0.0000001;
    /// opts.fill_gaps = true;
    /// opts.call_variants = true;
    /// opts.format = true;
    /// # let expected = kbo::MapOpts::default();
    /// # assert_eq!(opts.max_error_prob, expected.max_error_prob);
    /// # assert_eq!(opts.fill_gaps, expected.fill_gaps);
    /// # assert_eq!(opts.call_variants, expected.call_variants);
    /// # assert_eq!(opts.format, expected.format);
    /// # assert_eq!(opts.sbwt_build_opts.k, expected.sbwt_build_opts.k);
    /// # assert_eq!(opts.sbwt_build_opts.add_revcomp, expected.sbwt_build_opts.add_revcomp);
    /// # assert_eq!(opts.sbwt_build_opts.num_threads, expected.sbwt_build_opts.num_threads);
    /// # assert_eq!(opts.sbwt_build_opts.prefix_precalc, expected.sbwt_build_opts.prefix_precalc);
    /// # assert_eq!(opts.sbwt_build_opts.build_select, true);
    /// # assert_eq!(opts.sbwt_build_opts.mem_gb, expected.sbwt_build_opts.mem_gb);
    /// # assert_eq!(opts.sbwt_build_opts.dedup_batches, expected.sbwt_build_opts.dedup_batches);
    /// # assert_eq!(opts.sbwt_build_opts.temp_dir, expected.sbwt_build_opts.temp_dir);
    /// ```
    ///
    fn default() -> MapOpts {
        MapOpts {
            max_error_prob: 0.0000001,
            fill_gaps: true,
            call_variants: true,
            format: true,
            sbwt_build_opts: index::BuildOpts { build_select: true, ..Default::default() },
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

/// Calls variants between a query and a reference sequence.
///
/// Builds an SBWT index from `ref_seq` with parameters and resources from
/// `call_opts.sbwt_build_opts` and compares the index against `sbwt_query` and
/// `lcs_query` to determine substitutions, insertions, and deletions in the
/// query. This function also requires the `call_opts.max_error_prob` option to
/// [derandomize] the matching statistics vectors.
///
/// Returns a vector containing the identified [variants](variant_calling::Variant).
///
/// # Examples
/// ```rust
/// use kbo::call;
/// use kbo::build;
/// use kbo::index::BuildOpts;
/// use kbo::CallOpts;
/// use kbo::variant_calling::Variant;
///
/// let reference = b"TCGTGGATCGATACACGCTAGCAGGCTGACTCGATGGGATACTATGTGTTATAGCAATTCGGATCGATCGA";
/// let query =     b"TCGTGGATCGATACACGCTAGCCTGACTCGATGGGATACCATGTGTTATAGCAATTCCGGATCGATCGA";
///
/// let mut call_opts =  CallOpts::default();
/// call_opts.sbwt_build_opts.k = 20;
/// call_opts.max_error_prob = 0.001;
///
/// let (sbwt_query, lcs_query) = build(&[query.to_vec()], call_opts.sbwt_build_opts.clone());
///
/// let variants = call(&sbwt_query, &lcs_query, reference, call_opts);
///
/// // `variants` is a vector with elements:
/// // - Variant { query_pos: 22, query_chars: [65, 71, 71], ref_chars: [] }
/// // - Variant { query_pos: 42, query_chars: [84], ref_chars: [67] }
/// // - Variant { query_pos: 60, query_chars: [], ref_chars: [67] }
///
/// # assert_eq!(variants[0], Variant { query_pos: 22, query_chars: [65, 71, 71].to_vec(), ref_chars: [].to_vec() });
/// # assert_eq!(variants[1], Variant { query_pos: 42, query_chars: [84].to_vec(), ref_chars: [67].to_vec() } );
/// # assert_eq!(variants[2], Variant { query_pos: 60, query_chars: [].to_vec(), ref_chars: [67].to_vec() });
/// ```
///
pub fn call(
    sbwt_query: &SbwtIndexVariant,
    lcs_query: &sbwt::LcsArray,
    ref_seq: &[u8],
    call_opts: CallOpts,
) -> Vec<variant_calling::Variant> {
    let (sbwt_ref, lcs_ref) = index::build_sbwt_from_vecs(&[ref_seq.to_vec()], &Some(call_opts.sbwt_build_opts));

    match sbwt_query {
        SbwtIndexVariant::SubsetMatrix(ref sbwt_query) => {
            match sbwt_ref {
                SbwtIndexVariant::SubsetMatrix(ref sbwt_ref) => {
                    assert!(sbwt_ref.k() == sbwt_query.k());
                    assert!(sbwt_ref.alphabet() == sbwt_query.alphabet());
                    variant_calling::call_variants(
                        sbwt_query,
                        lcs_query,
                        sbwt_ref,
                        &lcs_ref,
                        ref_seq,
                        call_opts.max_error_prob,
                    )
                },
            }
        },
    }
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
/// use kbo::MatchOpts;
///
/// let reference: Vec<Vec<u8>> = vec![vec![b'A',b'A',b'A',b'G',b'A',b'A',b'C',b'C',b'A',b'-',b'T',b'C',b'A',b'G',b'G',b'G',b'C',b'G']];
/// let mut opts = BuildOpts::default();
/// opts.k = 3;
/// let (sbwt, lcs) = build(&reference, opts);
///
/// let query = vec![b'G',b'T',b'G',b'A',b'C',b'T',b'A',b'T',b'G',b'A',b'G',b'G',b'A',b'T'];
///
/// let ms_vectors = matches(&query, &sbwt, &lcs, MatchOpts::default());
/// // `ms_vectors` has ['-','-','-','-','-','-','-','-','-','M','M','M','-','-']
/// # assert_eq!(ms_vectors, vec!['-','-','-','-','-','-','-','-','-','M','M','M','-','-']);
/// ```
///
pub fn matches(
    query_seq: &[u8],
    sbwt: &SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
    match_opts: MatchOpts,
) -> Vec<char> {
    let (k, threshold) = match sbwt {
        SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
            (sbwt.k(), derandomize::random_match_threshold(sbwt.k(), sbwt.n_kmers(), 4_usize, match_opts.max_error_prob))
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
/// Run the full algorithm
/// ```rust
/// use kbo::build;
/// use kbo::map;
/// use kbo::index::BuildOpts;
/// use kbo::MapOpts;
///
/// let query: Vec<Vec<u8>> = vec![vec![b'A',b'A',b'A',b'G',b'A',b'A',b'C',b'C',b'A',b'-',b'T',b'C',b'A',b'G',b'G',b'G',b'C',b'G']];
/// let mut opts = BuildOpts::default();
/// opts.k = 3;
/// opts.build_select = true;
/// let (sbwt_query, lcs_query) = build(&query, opts.clone());
///
/// let mut map_opts = MapOpts::default();
/// map_opts.sbwt_build_opts = opts;
///
/// let reference = vec![b'G',b'T',b'G',b'A',b'C',b'T',b'A',b'T',b'G',b'A',b'G',b'G',b'A',b'T'];
///
/// let alignment = map(&reference, &sbwt_query, &lcs_query, map_opts);
/// // `ms_vectors` has [45,45,45,45,45,45,45,45,45,65,71,71,45,45]
/// # assert_eq!(alignment, vec![45,45,45,45,45,45,45,45,45,65,71,71,45,45]);
/// ```
///
/// Skip gap filling and variant calling, outputting '-'s instead.
/// ```rust
/// use kbo::build;
/// use kbo::map;
/// use kbo::index::BuildOpts;
/// use kbo::MapOpts;
///
/// let reference = b"CGTTGACTCTAGGTGCCTGGGTTCTCAGAGCTGGGC".to_vec();
/// let query =     b"CGTTGACTGGTGCCTGGGTTCTCAGAGCTGGGC".to_vec();
/// let expected =  b"CGTTGACT---GGTGCCTGGGTTCTCAGAGCTGGGC".to_vec();
///
/// let mut opts = BuildOpts::default();
/// opts.k = 7;
/// opts.build_select = true;
///
/// let mut map_opts = MapOpts::default();
/// map_opts.fill_gaps = false;
/// map_opts.call_variants = false;
/// map_opts.max_error_prob = 0.1;
/// map_opts.sbwt_build_opts = opts.clone();
///
/// let (sbwt_query, lcs_query) = build(&[query], opts);
///
/// let alignment = map(&reference, &sbwt_query, &lcs_query, map_opts);
/// // `alignment` has CGTTGACT---GGTGCCTGGGTTCTCAGAGCTGGGC
/// # assert_eq!(alignment, expected);
/// ```
///
/// Output the internal representation of [translate]
/// ```rust
/// use kbo::build;
/// use kbo::map;
/// use kbo::index::BuildOpts;
/// use kbo::MapOpts;
///
/// let reference = b"CGTTGACTCTAGGTGCCTGGGTTCTCAGAGCTGGGC".to_vec();
/// let query =     b"CGTTGACTGGTGCCTGGGTTCTCAGAGCTGGGC".to_vec();
/// let expected =  b"MMMMMMMM---MMMMMMMMMMMMMMMMMMMMMMMMM".to_vec();
///
/// let mut opts = BuildOpts::default();
/// opts.k = 7;
/// opts.build_select = true;
///
/// let mut map_opts = MapOpts::default();
/// map_opts.fill_gaps = false;
/// map_opts.call_variants = false;
/// map_opts.format = false;
/// map_opts.max_error_prob = 0.1;
/// map_opts.sbwt_build_opts = opts.clone();
///
/// let (sbwt_query, lcs_query) = build(&[query], opts);
///
/// let alignment = map(&reference, &sbwt_query, &lcs_query, map_opts);
/// // `alignment` has MMMMMMMM---MMMMMMMMMMMMMMMMMMMMMMMMM
/// # assert_eq!(alignment, expected);
/// ```
///
pub fn map(
    ref_seq: &[u8],
    query_sbwt: &SbwtIndexVariant,
    query_lcs: &sbwt::LcsArray,
    map_opts: MapOpts,
) -> Vec<u8> {
    let (k, threshold) = match query_sbwt {
        SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
            if map_opts.call_variants {
                assert!(sbwt.k() == map_opts.sbwt_build_opts.k);
            }
            (sbwt.k(), derandomize::random_match_threshold(sbwt.k(), sbwt.n_kmers(), 4_usize, map_opts.max_error_prob))
        },
    };

    let noisy_ms = index::query_sbwt(ref_seq, query_sbwt, query_lcs);
    let derand_ms = derandomize::derandomize_ms_vec(&noisy_ms.iter().map(|x| x.0).collect::<Vec<usize>>(), k, threshold);

    let translation = translate::translate_ms_vec(&derand_ms, k, threshold);

    let mut call_opts = CallOpts{ sbwt_build_opts: map_opts.sbwt_build_opts, ..Default::default() };
    call_opts.max_error_prob = map_opts.max_error_prob;

    let refined = if map_opts.fill_gaps {
        gap_filling::fill_gaps(&translation, &noisy_ms, ref_seq, query_sbwt, threshold, map_opts.max_error_prob)
    } else {
        translation
    };

    let with_variants = if map_opts.call_variants {
        let variants = call(query_sbwt, query_lcs, ref_seq, call_opts);
        translate::add_variants(&refined, &variants)
    } else {
        refined
    };

    if map_opts.format {
        format::relative_to_ref(ref_seq, &with_variants)
    } else {
        with_variants.iter().map(|x| *x as u8).collect()
    }
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
/// ```rust
/// use kbo::build;
/// use kbo::find;
/// use kbo::FindOpts;
/// use kbo::index::BuildOpts;
/// use kbo::format::RLE;
///
/// let gene1: Vec<u8> = b"ATGGCTGTTCCATCATCAAAAGAAGAGTTAATTAAAGCTATTAATAGTAATTTTTCTTTATTAAATAAGAAGCTAGAATCTATTACGCCCCAACTCGCCTTTGAACCTCTATTGGAAGGGCACGCGAAGGGGACTACGATTAGCGTAGCGAATCTGGTTTCCTATCTGATTGGCTGGGGAGAGCTGGTGTTACACTGGCATGACCAAGAGGCAAAAGGAAAAACTATTATTTTTCCTGAGGAAGGATTTAAATGGAATGAATTGGGGCGTTTAGCACAGAAATTCTACCGTGACTATGAGGATATTACAGAGTACGAAGTTTTATTGGCACGGTTAAAGGAAAATAAGCAGCAACTCGTGGCTTTGATTGAACGATTCAGTAACGACGAGCTTTACGGTAAACCTTGGTATAATAAATGGACCCGAGGTCGTATGATTCAATTTAATACCGCCTCGCCTTATAAAAATGCTTCGGGGAGGTTAAATAAACTGCAGAAATGTCTTGCAGAATAG".to_vec();
/// let gene2_rc: Vec<u8> = b"CTACCCTACTATTTCGAGTGATTCAATCGTCTGGTTCACATAACCTACCACCTGTTCAAAATGCTTATCGACAAAAAAATGATCGGCAGCAGGAAATATAATAGTCCGCGTCTTTCGTGTGGTGAATTTTTCCCATGCAAGTAATTCATCCTGCATTACCAGATTGTCAGCATCGCCATGAAATAGCACGATCGGACAGGTTAATGTGCGCGCCTTGGCCTGAAATACATACTGCTCATAGAGCCGATAATCGTTTTTAATGATGGGGGTGAAAATTGTCATTAACTCTTTATTACGAAAGACATCAACCGGAGTTCCGCCCAGCTTGACGATCTCTTCCATAAACGCCTGATCGGGCAAGGTATGCAGTATTACTTCATGAGAGGCCCGATCGGGTGGGCGACAGCCGGAAAAAAACAGCGCGCATGGCATGTCATGTCCATGATCGAGAATATAATGCACCAGTTCGAAGGCCATGATCCCTCCGAGACTATGCCCAAAAATGGCGTAGTCTCCACCTGTGTAGTGTTTCACAAATTGTTGATAAAGGTCAGCGACGGCATCCACCATCGTAAGACACAGCGGCTGGCGTATTCTAGTTCCCCTCCCCGCAGGTTCTAAAGGCCGCAAAGTAATATTGTCCGACAGCACGCTACGCCATTTATAATACATGGCGGCAGAACCACCTGAATATGGCAAACAATACAAACTGATATTACTCAT".to_vec();
/// let reference: Vec<Vec<u8>> = vec![gene1, gene2_rc];
///
/// let query: Vec<u8> = b"ATGGCTGTTCCATCATCAAAAGAAGAGTTAATTAAAGCTATTAATAGTAATTTTTCTTTATTAAATAAGAAGCTAGACTCTATTACGCCCCAACTCGCCTTTGAACCTCTATTGGAAGGGCACGCGAAGGGGACTACGATTAGCGTAGCGAATCTGGTTTCCTATCTGATTGGCTGGGGAGAGCTGGTGTTACACTGGCATGACCAAGAGGCAAAAGGAAAAACTATTATTTTTCCTGAGGAAGGATTTAAATGGAATGAATTGGGGCGTTTAGCACAGAAATTCTACCGTGACTATGAGGATATTACAGAGTACGAAGTTTTATTGGCACGGTTAAAGGAAAATAAGCAGCAACTCGTGGCTTTGATTGAACGATTCAGTAACGACGAGCTTTACGGTAAACCTTGGTATAATAAATGGACCCGAGGTCGTATGATTCAATTTAATACCGCCTCGCCTTATAAAAATGCTTCGGGGAGGTTAAATAAACTGCAGAAATGTCTTGCAGAATAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCTACCCTACTATTTCGAGTGATTCAATCGTCTGGTTCACATAACCTACCACCTGTTCAAAATGCTTATCGACAAAAAAATGATCGGCAGCAGGAAATATAATAGTCCGCGTCTTTCGTGTGGTGAATTTTTCCCATGCAAGTAATTCATCCTGCATTACCAGATTGTCAGCATCGCCATGAAATAGCACGATCGGACAGGTTAATGTGCGCGCCTTGGCCTGAAATACATACTGCTCATAGAGCCGATAATCGTGTGTAATGATGGGGGTGAAAATTGTCATTAACTCTTTATTACGAAAGACATCAACCGGAGTTCCGCCCAGCTTGACGATCTCTTCCATAAACGCCTGATCGGGCAAGGTATGCATTTTTTTTTTTTTTTTTTTTTTTTTTGTATTACTTCATGAGAGGCCCGATCGGGTGGGCGACAGCCGGAAAAAAACAGCGCGCATGGCATGTCATGTCCATGATCGAGAATATAATGCACCAGTTCGAAGGCCATGATCCCTCCGAGACTATGCCCAAAAATGGCGTAGTCTCCACCTGTGTAGTGTTTCACAAATTGTTGATAAAGGTCAGCGACGGCATCCACCATCGTAAGACACAGCGGCTGGCGTATTCTAGTTCCCCTCCCCGCAGGTTCTAAAGGCCGCAAAATAATATTGCGACAGCACGCTACGCCATTTATAATACATGGCGGCAGAACCACCTGAATATGGCAAACAATACAAACTGATATTACTCAT".to_vec();
///
/// let mut opts = BuildOpts::default();
/// opts.k = 31;
/// let (sbwt, lcs) = build(&reference, opts);
///
/// let mut find_opts = FindOpts::default();
/// find_opts.max_gap_len = 50;
///
/// let local_alignments = find(&query, &sbwt, &lcs, find_opts);
///
/// // `local_alignments` has:
/// // - RLE { start: 0, end: 513, matches: 512, mismatches: 1, jumps: 0, gap_bases: 0, gap_opens: 0 },
/// // - RLE { start: 593, end: 1340, matches: 709, mismatches: 0, jumps: 0, gap_bases: 38, gap_opens: 3 }
///
/// # assert_eq!(local_alignments, [RLE { start: 0, end: 513, matches: 512, mismatches: 1, jumps: 0, gap_bases: 0, gap_opens: 0 }, RLE { start: 593, end: 1340, matches: 709, mismatches: 0, jumps: 0, gap_bases: 38, gap_opens: 3 }]);
/// ```
///
pub fn find(
    query_seq: &[u8],
    sbwt: &SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
    find_opts: FindOpts,
) -> Vec<format::RLE> {
    let match_opts = MatchOpts { max_error_prob: find_opts.max_error_prob };
    let aln = matches(query_seq, sbwt, lcs, match_opts);
    if find_opts.max_gap_len > 0 {
        format::run_lengths_gapped(&aln, find_opts.max_gap_len)
    } else {
        format::run_lengths(&aln)
    }
}
