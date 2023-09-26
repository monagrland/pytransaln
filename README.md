Translation-guided alignment of nucleotide sequences
====================================================

Existing tools for translation-guided codon alignment may no longer be
accessible, e.g. [TranslatorX](https://doi.org/10.1093/nar/gkq291) or
[pal2nal](https://www.bork.embl.de/pal2nal/), or may need to be ported to new
language or dependency versions, e.g.
[TransAlign](https://uol.de/systematik-evolutionsbiologie/programme).

This is a reimplementation of some features of the above programs to perform
simple translation-guided nucleotide (codon) alignments.

Reading frame can be manually specified or guessed with a heuristic. Genetic
code must be manually specified; heuristic to guess genetic code is not
implemented.

Frameshifted sequences are treated here with a heuristic based on them having
more stop codons. For a more careful alignment, or for sets with many
frameshifted sequences, use [MACSE](https://www.agap-ge2pop.org/macse/)
instead, however MACSE is quite slow for de novo alignments and is probably
overkill for most "normal" datasets without many frameshifts or pseudogenes.


## Installation

Use pip to install from the source folder; recommended to install into a
virtualenv.

```bash
pip install .
```


## Usage

Requires [MAFFT](https://mafft.cbrc.jp/alignment/software/) >=6.811; tested
with v7.520.

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
