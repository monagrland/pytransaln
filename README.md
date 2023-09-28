Translation-guided alignment of nucleotide sequences
====================================================

Existing tools for translation-guided codon alignment may no longer be
accessible, e.g. [TranslatorX](https://doi.org/10.1093/nar/gkq291) or
[pal2nal](https://www.bork.embl.de/pal2nal/), or may need to be ported to new
language or dependency versions, e.g.
[TransAlign](https://uol.de/systematik-evolutionsbiologie/programme).

This is a reimplementation of some features of the above programs to perform
simple translation-guided nucleotide (codon) alignments, and to screen for
pseudogenes with frameshift indels or non-sense substitutions.

The tool can be used to perform alignment or simply report sequence statistics
and flag potential pseudogenes.


## How the alignment works

Reading frame can be manually specified or guessed with a heuristic. Genetic
code must be manually specified; heuristic to guess genetic code is not yet
implemented.

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
* Translated 'good' sequences are aligned with MAFFT; amino acid translation is
  then back-translated to a codon alignment
* Nucleotide sequences of putative pseudogenes are then aligned against the
  reference 'good' alignment with MAFFT `--add` option
* Likely frameshift positions in putative pseudogenes are reported from the
  positional map of the reference-guided alignment


## Assumptions

* Input sequences are homologous
* Input sequences are protein coding sequences without introns or untranslated regions
* Input sequences are long enough that wrong reading frame will be evident in excessive stop codons
* If pseudogenes are present, majority of sequences are not pseudogenes
* Sequences all use the same genetic code

For a more careful alignment, or for sequence sets with many frameshifted
sequences, use [MACSE](https://www.agap-ge2pop.org/macse/) instead, however
MACSE is quite slow for de novo alignments and is probably overkill for most
"normal" datasets without many frameshifts or pseudogenes.


## Installation

Use pip to install from the source folder; recommended to install into a
virtualenv.

```bash
pip install .
```

External dependencies not installed via pip, should be in path:
* [MAFFT](https://mafft.cbrc.jp/alignment/software/) >=6.811; tested with v7.520.


## Usage

See help message for details

```bash
pytransaln --help
```

Recommended to inspect alignment afterwards or apply quality checks with other
programs such as [trimAl](http://trimal.cgenomics.org/).

To view alignments on the command line you can use
[alv](https://github.com/arvestad/alv) and pipe to less with the `-R` option:

```bash
alv -t codon -l alignment.fasta | less -R
```


## Future enhancements

In order of priority

- [x] Diagnostic plots: stop codons per frame; min stop codons per seq
- [ ] Screen sequences with HMM profile of protein sequence
- [ ] Identify sequences with wrong frame or frameshifts
- [ ] User-supplied input amino acid alignment
- [x] Identify likely frameshift positions from MAFFT .map file
- [ ] Add pre and post frame sequence to alignment
- [ ] Guess genetic code
- [ ] Translate 6 frames
