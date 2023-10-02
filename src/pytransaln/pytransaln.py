#!/usr/bin/env python3

from pytransaln.stats import stats
from pytransaln.align import align
from pytransaln import __version__

import argparse
import logging
import sys

logging.basicConfig(
    format="%(asctime)s : %(levelname)s : %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="Translation-guided nucleotide alignment and reading frame statistics of coding sequences"
    )
    parser.add_argument(
        "--input",
        help="Path to input file with unaligned nucleotide sequences, Fasta format",
    )
    parser.add_argument(
        "--code",
        default=5,
        type=int,
        help="Genetic code to use for all sequences, NCBI translation table number (except stopless codes 27, 28, 31)",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=__version__,
    )
    subparsers = parser.add_subparsers(required=True)
    parser_align = subparsers.add_parser(
        "align", help="Translation-guided nucleotide (codon) alignment"
    )
    parser_align.set_defaults(func=align)
    parser_stats = subparsers.add_parser(
        "stats", help="Report summary stats on reading frames only"
    )
    parser_stats.set_defaults(func=stats)
    # Alignment
    parser_align.add_argument(
        "--how",
        default="u",
        help="""
        How to choose reading frame: 'each' - find reading frame that minimizes
        stop codons for each individual sequence; may result in more than one
        possible frame per sequence; 'cons' - find frame that minimizes stop
        codons across all sequences and apply that frame too all sequences;
        'user' - user specified reading frame at option --frame
        """,
    ),
    parser_align.add_argument(
        "--maxstops",
        default=0,
        type=int,
        help="Max stop codons to allow in 'good' alignment; nt sequences over this threshold in all frames will be written to --out_bad",
    )
    parser_align.add_argument(
        "--frame",
        default=0,
        type=int,
        help="Reading frame offset to apply to all sequences, must be 0, 1, or 2; overridden by --how each or --how cons",
    )
    parser_align.add_argument(
        "--aligner",
        default="mafft",
        help="Alignment program to use (only MAFFT implemented at the moment)",
    )
    parser_align.add_argument(
        "--out_aa",
        default="test.aa.fasta",
        help="Path to write translated AA sequences",
    )
    parser_align.add_argument(
        "--out_bad",
        default="test.bad.nt.fasta",
        help="Path to write sequences with too many stop codons",
    )
    parser_align.add_argument(
        "--out_aln_aa",
        default="test.aln.aa.fasta",
        help="Path to write aligned amino acid sequences",
    )
    parser_align.add_argument(
        "--out_aln_nt",
        default="test.aln.nt.fasta",
        help="Path to write aligned nucleotide sequences",
    )
    parser_align.add_argument(
        "--out_aln_nt_aug",
        default="test.aln.nt.aug.fasta",
        help="Path to write aligned nucleotide sequences with likely frameshifted sequences added",
    )
    parser_align.add_argument(
        "--out_bad_fs_report",
        default="test.bad.nt.frameshifts.tsv",
        help="Path to write report on likely frameshifts in 'bad' nucleotide sequences",
    )
    parser_align.add_argument(
        "--threads",
        default=1,
        type=int,
        help="Number of threads to pass to alignment program",
    )
    # Stats
    parser_stats.add_argument(
        "--hmm",
        help="Optional HMM to screen translated sequences; only first model will be read from file",
    )
    parser_stats.add_argument(
        "--out_hmmsearch",
        default="test.hmmsearch_tblout.txt",
        help="Path to write tabular output from hmmsearch",
    )
    parser_stats.add_argument(
        "--out_hist_hmm",
        default="test.hist_hmm_scores.png",
        help="Path to plot histogram of HMM bit scores",
    )
    parser_stats.add_argument(
        "--out_screened",
        default="test.screened_ok.fasta",
        help="Path to write sequences that passed screening, Fasta format",
    )
    parser_stats.add_argument(
        "--out_stats",
        default="test.stopcodon_stats.tsv",
        help="Path to write per-frame stop codon statistics",
    )
    parser_stats.add_argument(
        "--out_hist_spf",
        default="test.hist_stops_perframe.png",
        help="Path to plot histogram of stops per reading frame",
    )
    parser_stats.add_argument(
        "--out_hist_mins",
        default="test.hist_minstops_perseq.png",
        help="Path to plot histogram of minimum stop codons per sequence",
    )
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
