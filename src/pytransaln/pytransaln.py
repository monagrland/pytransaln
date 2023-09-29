#!/usr/bin/env python3

from pytransaln.frameshifts import report_frameshifts
from pytransaln.translate import translate_3_frames, translate_1_frame, onebestframe, guessframe
import argparse
import logging
import sys
import pyhmmer
import pandas as pd
import matplotlib.pyplot as plt

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


def summarize_framestats(trseq):
    """Tabulate number stop codons per frame from translate_3_frames output"""
    summary = [
        {"seq_id": i, "frame": frame, "stops": trseq[i][frame].seq.count("*")}
        for i in trseq
        for frame in trseq[i]
    ]
    return pd.DataFrame.from_dict(summary)


def seqrecords2sequenceblock(seqrecords, alphabet=pyhmmer.easel.Alphabet.amino()):
    """Convert list of Biopython SeqRecords to Easel Sequenceblock

    Parameters
    ----------
    seqrecords : list
        The .id attribute of each record will be used to populate .name
        attribute of the sequences in Sequenceblock.
    alphabet : pyhmmer.easel.Alphabet object

    Returns
    -------
    pyhmmer.easel.TextSequenceBlock
    """
    seqblock = pyhmmer.easel.TextSequenceBlock(
        [
            pyhmmer.easel.TextSequence(sequence=str(i.seq), name=i.id.encode())
            for i in seqrecords
        ]
    )
    seqblock = seqblock.digitize(alphabet)
    return seqblock


def summarize_framestats_with_hmm(trseq, hmmfile, outfile=None):
    """Tabulate stop codons and HMM score of three-frame translation

    Parameters
    ----------
    trseq : dict
        Output of translate_3_frames
    hmmfile : str
        Path to HMM file; only the first model in file will be used.
    outfile : str
        Path to write HMM results in tblout format

    Returns
    -------
    pd.DataFrame
    """
    seqlist = [trseq[i][frame] for i in trseq for frame in trseq[i]]
    seqblock = seqrecords2sequenceblock(
        seqlist, alphabet=pyhmmer.easel.Alphabet.amino()
    )
    with pyhmmer.plan7.HMMFile(hmmfile) as hmm_file:
        hmm = hmm_file.read()
    pipeline = pyhmmer.plan7.Pipeline(hmm.alphabet)
    hits = pipeline.search_hmm(hmm, seqblock)
    if outfile:
        with open(outfile, "wb") as fh:
            hits.write(fh, format="targets")
    id2score = {i.name.decode(): i.score for i in hits}
    summary = [
        {
            "seq_id": i,
            "frame": frame,
            "stops": trseq[i][frame].seq.count("*"),
            "hmm_score": id2score[trseq[i][frame].id]
            if trseq[i][frame].id in id2score
            else None,
        }
        for i in trseq
        for frame in trseq[i]
    ]
    return pd.DataFrame.from_dict(summary)


def hist_stops_per_frame(df):
    """Plot histogram of the number of stop codons facetted by reading frame

    If sequences are amplified by the same PCR primers, we expect them all to
    be in the same frame. Most sequences in the correct reading frame should
    have zero stop codons.

    Possible exceptions: Wrong genetic code used; amplified sequence
    encompasses introns or untranslated regions; sequences mostly pseudogenes
    or non-coding.

    Parameters
    ----------
    df : pandas.DataFrame
        Output from summarize_framestats()

    Returns
    -------
    fig, axs
    """
    fig, axs = plt.subplots(3, 1, sharex=True, sharey=True, layout="constrained")
    breaks = range(0, df["stops"].max() + 1)
    for frame in [0, 1, 2]:
        axs[frame].hist(
            x=df.query(f"frame == {frame}")["stops"],
            bins=breaks,
        )
        axs[frame].set_title(f"Frame offset {str(frame)}")
        axs[frame].set_ylabel("Count")
    axs[2].set_xlabel("Stop codons per sequence")
    axs[2].set_xticks(breaks)
    return fig, axs


def hist_minstops_per_seq(df):
    """Plot histogram of minimum number of stops per sequence

    Parameters
    ----------
    df : pandas.DataFrame
        Output from summarize_framestats()

    Returns
    -------
    fig, axs
    """
    minstops = df.groupby("seq_id")[["stops"]].min()
    fig, axs = plt.subplots(1)
    breaks = range(0, minstops["stops"].max() + 1)
    axs.hist(minstops["stops"], bins=breaks)
    axs.set_ylabel("Count")
    axs.set_xlabel("Minimum number of stop codons")
    axs.set_xticks(breaks)
    return fig, axs


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


def stats(args):
    nt = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
    logger.info("%d nucleotide sequences in input", len(nt))
    seq2code = {i: args.code for i in nt}
    trseq = translate_3_frames(nt, seq2code)
    if args.hmm:
        logger.info("Using HMM model in %s to screen translations", args.hmm)
        df = summarize_framestats_with_hmm(trseq, args.hmm, args.out_hmmsearch)
    else:
        df = summarize_framestats(trseq)
    logger.info("Writing summary stats to %s", args.out_stats)
    df.to_csv(args.out_stats, sep="\t", index=False)
    # Histograms
    logger.info("Plotting histograms to %s and %s", args.out_hist_spf, args.out_hist_mins)
    hist_spf_fig, hist_spf_axs = hist_stops_per_frame(df)
    hist_spf_fig.savefig(args.out_hist_spf)
    hist_mins_fig, hist_mins_axs = hist_minstops_per_seq(df)
    hist_mins_fig.savefig(args.out_hist_mins)


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
