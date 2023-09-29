#!/usr/bin/env python3

from pytransaln.frameshifts import report_frameshifts
from pytransaln.translate import translate_3_frames, translate_1_frame, onebestframe, guessframe
from pytransaln.stats import stats

import argparse
import logging
import sys
import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO
from subprocess import run
from collections import defaultdict

logging.basicConfig(
    format="%(asctime)s : %(levelname)s : %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)


def trdict2seqlist(trdict):
    out = [trdict[i][frame] for i in trdict for frame in trdict[i]]
    aa2nt = {trdict[i][frame].id: i for i in trdict for frame in trdict[i]}
    return out, aa2nt


def yield_codons(seq):
    """Generator to yield codon triplets from nucelotide sequence

    Parameters
    ----------
    seq : Bio.SeqRecord.SeqRecord or str
        Input nucleotide sequence

    Returns
    -------
    Bio.SeqRecord.SeqRecord, Seq, or str
        Codon triplets
    """
    for i in range(int(len(seq) / 3)):
        yield seq[i * 3 : i * 3 + 3]
    if len(seq) % 3 != 0:  # trailing untranslated sequence
        yield seq[int(len(seq) / 3) * 3 :]


def aa_aln_to_nt_aln(aa, nt, frame=0):
    """Align nucleotide sequence against aligned amino acid sequence

    Does not check for correctness of the translation. Unlike Bio.codonalign,
    reading frame offset is taken into account.

    Parameters
    ----------
    aa : Bio.SeqRecord.SeqRecord
        Aligned amino acid sequence, gaps as '-'
    nt : Bio.SeqRecord.SeqRecord
        Unaligned nucleotide sequence, assumed to have no gaps!
    frame : int
        Reading frame offset for the nucleotide sequence; 0, 1, or 2

    Returns
    -------
    (str, str, str)
        Tuple of nucleotide sequences representing the initial offset base(s),
        the nucleotide sequence aligned in codon blocks, and any trailing
        unaligned bases
    """
    yc = yield_codons(nt.seq[frame:])
    pre = str(nt.seq[0:frame])
    out = ""
    for i in aa.seq:
        if i == "-":
            out += "---"
        else:
            out += str(next(yc))
    if len(nt.seq[frame:]) % 3 != 0:
        post = str(next(yc))
    else:
        post = ""
    return pre, out, post


def align(args):
    if args.frame not in [0, 1, 2]:
        raise ValueError("Frame must be 0, 1, or 2")

    if args.code in [27, 28, 31]:
        raise ValueError("Please choose a non-ambiguous genetic code")

    nt = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
    logger.info(f"{str(len(nt))} nucleotide sequences to align")

    too_many_stops = None
    if args.how.startswith("e"):
        logger.info(
            "Guessing reading frame for each sequence by minimizing stop codons"
        )
        seq2code = {i: args.code for i in nt}
        tr, too_many_stops = guessframe(
            seqdict=nt, codes=seq2code, maxstops=args.maxstops
        )
        seq2frame = {i: list(tr[i].keys())[0] for i in tr}
    elif args.how.startswith("c"):
        logger.info(
            "Find one reading frame for all sequence that minimizes total stop codons"
        )
        seq2code = {i: args.code for i in nt}
        tr, too_many_stops = onebestframe(
            seqdict=nt, codes=seq2code, maxstops=args.maxstops
        )
        seq2frame = {i: list(tr[i].keys())[0] for i in tr}
    else:
        logger.info(
            f"Applying same reading frame offset {str(args.frame)} for all sequences"
        )
        # Single reading frame for all sequences
        seq2frame = {i: args.frame for i in nt}
        seq2code = {i: args.code for i in nt}
        tr, too_many_stops = translate_1_frame(
            seqdict=nt, frames=seq2frame, codes=seq2code, maxstops=args.maxstops
        )

    if too_many_stops:
        logger.info(
            f"{str(len(too_many_stops))} sequences with > {str(args.maxstops)} stop codons"
        )
        if len(too_many_stops) >= 0.5 * (len(nt)):
            logger.info(
                "More than 50% of sequences have too many stop codons; check genetic code and sequence orientation?"
            )
    logger.info(f"{str(len(tr))} sequences for initial alignment")

    aa, aa2nt = trdict2seqlist(tr)

    with open(args.out_aa, "w") as fh:
        SeqIO.write(aa, fh, "fasta")

    # read aa alignment
    logger.info("Aligning with MAFFT")
    cmd = ["mafft", "--thread", str(args.threads), args.out_aa]
    logger.info("Command: %s", " ".join(cmd))
    mafft_job = run(cmd, capture_output=True)
    logger.debug(mafft_job.stderr.decode())
    traln = SeqIO.to_dict(SeqIO.parse(StringIO(mafft_job.stdout.decode()), "fasta"))
    with open(args.out_aln_aa, "w") as fh:
        SeqIO.write(list(traln.values()), fh, "fasta")

    # align nt to aa
    ntaln = []
    for i in traln:
        pre, mid, post = aa_aln_to_nt_aln(traln[i], nt[aa2nt[i]], seq2frame[aa2nt[i]])
        ntaln.append(SeqRecord(Seq(mid), id=aa2nt[i], name=aa2nt[i]))
    with open(args.out_aln_nt, "w") as fh:
        SeqIO.write(ntaln, fh, "fasta")

    # add nt sequences with too many stop codons to the "clean" alignment
    if too_many_stops:
        logger.info("Adding putative pseudogenes to initial alignment")
        with open(args.out_bad, "w") as fh:
            SeqIO.write([nt[i] for i in too_many_stops], fh, "fasta")
        cmd = [
            "mafft",
            "--add",
            args.out_bad,
            "--mapout",
            "--thread",
            str(args.threads),
            args.out_aln_nt,
        ]
        logger.info("Command: %s", " ".join(cmd))
        mafft_add = run(cmd, capture_output=True)
        logger.debug(mafft_add.stderr.decode())
        with open(args.out_aln_nt_aug, "w") as fh:
            fh.write(mafft_add.stdout.decode())
        mapout = args.out_bad + ".map"
        frameshifts = report_frameshifts(mapout)
        for i in frameshifts:
            logger.info(
                "Sequence %s has %d likely frameshifts", i, len(frameshifts[i])
            )
        dfs = {
            i: pd.DataFrame(frameshifts[i]) for i in frameshifts
        }  # TODO: coordinate columns as integers
        for i in dfs:
            dfs[i]["seq_id"] = i
        df = pd.concat(list(dfs.values()))
        df.to_csv(args.out_bad_fs_report, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", help="Input unaligned nucleotide sequences, Fasta format"
    )
    parser.add_argument(
        "--code",
        default=5,
        type=int,
        help="Genetic code to use for all sequences, NCBI translation table number (except stopless codes 27, 28, 31)",
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
        "--out_stats",
        default="test.stopcodon_stats.tsv",
        help="Path to write per-frame stop codon statistics",
    )
    parser_stats.add_argument(
        "--out_hist_spf",
        default="test.stopsperframe.png",
        help="Path to plot histogram of stops per reading frame",
    )
    parser_stats.add_argument(
        "--out_hist_mins",
        default="test.minstopsperseq.png",
        help="Path to plot histogram of minimum stop codons per sequence",
    )
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
