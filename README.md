Translation-guided alignment of nucleotide sequences
====================================================

Several established tools for translation-guided codon alignment are no longer
maintained or available for download, e.g.
[TranslatorX](https://doi.org/10.1093/nar/gkq291) or
[pal2nal](https://www.bork.embl.de/pal2nal/), or may need to be ported to new
language or dependency versions to work properly, e.g.
[transAlign](https://uol.de/systematik-evolutionsbiologie/programme).

This package reimplements some features of the above programs to perform simple
translation-guided nucleotide (codon) alignments, and to screen for pseudogenes
with frameshift indels or non-sense substitutions.

The tool can be used to perform alignment or simply report sequence statistics
and flag potential pseudogenes. The intended use case is to screen and align
collections of PCR-amplified coding sequences used for metabarcoding, e.g. the
mitochondrial cytochrome c oxidase subunit I (mtCOI) gene fragment.


## How the alignment works

Reading frame can be manually specified or guessed with a heuristic. Genetic
code must be manually specified; heuristic to guess genetic code is not yet
implemented.

Choose reading frames for translation:
* Reading frame can be chosen in one of three ways (specified to option `--how`):
  * User-defined frame offset applied to all sequences (`--how user`)
  * Apply same frame to all sequences, choose consensus frame that minimizes
    total number of stop codons across all sequences (`--how cons`)
  * Choose frame individually for each sequence that minimizes stop codons for
    that sequence; may result in ties where a sequence may have more than one
    'best' reading frame (`--how each`)
* Sequences that have more than the maximum allowed number of stop codons in
  any reading frame are flagged as putative pseudogenes. 
* The 'good' sequences are translated in the reading frame as chosen above.
* If there is more than one reading frame with zero stop codons, the two (or
  three) alternative translations are each pairwise aligned to the remaining
  sequences with an unambiguous best reading frame. The frame that has the
  highest total alignment score is chosen.
* Optional: If an HMM representing the target protein sequence is provided
  (option `--hmm`), the 'good' sequences will be screened against this HMM;
  sequences with outlier bit scores will not be used for the initial aligment
* Translated 'good' sequences are aligned with MAFFT; nucleotide sequences are
  then aligned as codons using this amino acid alignment

Dealing with pseudogenes/frameshifted sequences (adapted from transAlign, see
[Bininda-Edmonds, 2005](https://doi.org/10.1186/1471-2105-6-156)):
* Nucleotide sequences of putative pseudogenes are then aligned against the
  reference 'good' alignment with MAFFT `--add` option
* Likely frameshift positions in putative pseudogenes are reported from the
  positional map of the reference-guided alignment


## Assumptions

* Input sequences are homologous
* Input sequences are protein coding sequences without introns or untranslated
  regions
* Input sequences are long enough that wrong reading frame will be evident in
  excessive stop codons (warning if average sequence length is under 50)
* If pseudogenes are present, majority of sequences are not pseudogenes
  (warning if more than half of sequences have excessive stop codons)
* Sequences all use the same genetic code

For a more careful alignment, or for sequence sets with many frameshifted
sequences, use [MACSE](https://www.agap-ge2pop.org/macse/) instead, however
MACSE is quite slow for de novo alignments and is probably overkill for most
"normal" datasets where most sequences do not have frameshifts.


## Installation

Install from PyPI with pip, preferably into a virtualenv or Conda environment:

```bash
pip install pytransaln
```

External dependencies are not installed via pip, but should also be in path:
* [MAFFT](https://mafft.cbrc.jp/alignment/software/) >=6.811; tested with v7.520.


## Usage

See help message for details

```bash
pytransaln --help
```

It is recommended to inspect the alignment afterwards or apply quality checks
with other programs such as [trimAl](http://trimal.cgenomics.org/).

To view alignments on the command line you can use
[alv](https://github.com/arvestad/alv) and pipe to less with the `-R` option:

```bash
alv -t codon -l alignment.fasta | less -R
```


## Output alignment

Note that:
* Leading and trailing bases not contained in complete codons may be omitted
* Portions of the putative pseudogene sequences that were aligned in the second
  step to the initial 'good' alignment may be omitted from the final
  'augmented' alignment to preserve the reading frame. See the frameshift
  report file for details.


## Testing and benchmarking

Commands to run tests with example data (from the benchmark data sets
distributed with
[transAlign](https://uol.de/systematik-evolutionsbiologie/programme)) are in
the Makefile:

```bash
make help # list available commands
make benchmark # download test data and run alignments
make clean # delete benchmark output
```


## Future enhancements

In order of priority

- [ ] Add pre and post frame sequence back to alignment
- [ ] Guess genetic code
- [ ] Translate 6 frames
- [ ] User-supplied input amino acid alignment
