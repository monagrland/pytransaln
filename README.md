Translation-guided alignment of nucleotide sequences
====================================================

Existing tools for translation-guided codon alignment may no longer be
accessible, e.g. [TranslatorX](https://doi.org/10.1093/nar/gkq291), or may need
to be ported to new language or dependency versions, e.g.
[TransAlign](https://uol.de/systematik-evolutionsbiologie/programme)

This is a reimplementation of some features of the above programs to perform
simple translation-guided alignments.

Reading frame can be manually specified or guessed with a heuristic. Genetic
code must be manually specified; heuristic to guess genetic code is not
implemented.

Frameshifted sequences are treated here with a heuristic based on them having
more stop codons. For a more careful alignment, or for sets with many
frameshifted sequences, use [MACSE](https://www.agap-ge2pop.org/macse/)
instead.
