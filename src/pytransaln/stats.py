#!/usr/bin/env python3

from pytransaln.translate import translate_3_frames
from pytransaln import __version__

import pyhmmer
import logging
import json
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from numpy import histogram

from Bio import SeqIO

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


def summarize_framestats_with_hmm(trseq, hmmfile, outfile=None, iqr_mult=1.5):
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
    df = pd.DataFrame.from_dict(summary)
    # Outlier HMM bit scores
    logger.info("Using outlier threshold multiplier of %.2f", iqr_mult)
    q1, q3 = df["hmm_score"].quantile([0.25, 0.75])  # NaN values are ignored
    iqr = q3 - q1
    ulim = q3 + iqr_mult * iqr
    llim = q1 - iqr_mult * iqr
    logger.info("Outlier thresholds for HMM bit scores: %d , %d", llim, ulim)
    df["hmm_ok"] = df["hmm_score"].apply(lambda x: x > llim and x < ulim)
    df["ok"] = df["hmm_ok"] & (df["stops"] == 0)
    return df, llim, ulim


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


def mqc_hist_stops_per_frame(df):
    """Dict of the number of stop codons facetted by reading frame for MultiQC

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
    dict
    """
    out = {
        "id" : "pytransaln_stops_per_frame",
        "section_name" : "Stop codons per sequence per reading frame",
        "description" : f"Counts of stop codons per sequence in each reading frame; if sequences are amplified by same PCR primers, expect them all to be in the same frame. Majority of sequences in that frame should have zero stop codons. Source: pytransaln v{__version__}",
        "plot_type" : "linegraph",
    }
    out['data'] = { frame : Counter(df.query(f"frame == {frame}")["stops"]) for frame in [0,1,2] }
    return out


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


def mqc_hist_minstops_per_seq(df):
    """Dict output of minimum number stops per sequence for MultiQC plot

    Parameters
    ----------
    df : pandas.DataFrame
        Output from summarize_framestats()

    Returns
    -------
    dict
    """
    out = {
        "id" : "pytransaln_minstops_per_seq",
        "section_name" : "Minimum stop codons per sequence",
        "description" : f"Number of stop codons in the reading frame with fewest stop codons, for each sequence. The majority should have 0 stops, if not, possible reasons include: wrong genetic code specified, sequence is reversed, sequence includes non-coding regions (e.g. introns). Source: pytransaln v{__version__}",
        "plot_type" : "linegraph",
    }
    minstops = df.groupby("seq_id")[["stops"]].min()
    minstops = dict(Counter(minstops['stops']))
    minstops = { int(i) : minstops[i] for i in minstops }
    out['data'] = { 'seqs' : minstops }
    return out


def hist_hmm_scores(df, vlines):
    """Plot histogram of minimum number of stops per sequence

    Parameters
    ----------
    df : pandas.DataFrame
        Output from summarize_framestats()
    vlines : list
        List of coordinates to draw vertical lines overlay to annotate plot

    Returns
    -------
    fig, axs
    """
    fig, ax = plt.subplots(1)
    histvals = ax.hist(df["hmm_score"])
    ax.vlines(vlines, ymin=0, ymax=histvals[0].max(), color="grey", linestyles="dashed")
    ax.set_ylabel("Count")
    ax.set_xlabel("HMM alignment bit score")
    return fig, ax


def mqc_hist_hmm_scores(df, vlines):
    """Dict output histogram of minimum number of stops per sequence for MultiQC

    Parameters
    ----------
    df : pandas.DataFrame
        Output from summarize_framestats()
    vlines : list
        List of coordinates to draw vertical lines overlay to annotate plot

    Returns
    -------
    dict
    """
    hist, bin_edges = histogram(df["hmm_score"].dropna(), bins="auto")
    out = {
        "id" : "pytransaln_hmm_scores_perseq",
        "section_name" : "HMM bit scores distribution",
        "description" : f"Bit scores of HMM alignments vs. translations, used to filter outlier sequences. Counts have been binned into histogram. Source: pytransaln v{__version__}",
        "plot_type" : "linegraph",
    }
    out['data'] = { 'seqs' : dict(zip(bin_edges.tolist(), hist.tolist())) } # serialize with .tolist() to change numpy int64 to python scalars, else json cannot dump
    return out


def stats(args):
    nt = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
    logger.info("%d nucleotide sequences in input", len(nt))
    seq2code = {i: args.code for i in nt}
    trseq = translate_3_frames(nt, seq2code, args.ignore_terminal_stop)
    if args.hmm:
        logger.info("Using HMM model in %s to screen translations", args.hmm)
        # TODO check for sequences which pass in more than one frame
        df, llim, ulim = summarize_framestats_with_hmm(
            trseq, args.hmm, args.out_hmmsearch, 1.5
        )
    else:
        df = summarize_framestats(trseq)
        df["ok"] = df["stops"] == 0
    logger.info("Writing summary stats to %s", args.out_stats)
    df.to_csv(args.out_stats, sep="\t", index=False)
    # Histograms
    logger.info(
        "Plotting histograms to %s and %s", args.out_hist_spf, args.out_hist_mins
    )
    hist_spf_fig, hist_spf_axs = hist_stops_per_frame(df)
    hist_spf_fig.savefig(args.out_hist_spf)
    dict_spf = mqc_hist_stops_per_frame(df)
    with open(args.out_mqc_spf, "w") as fh:
        json.dump(dict_spf, fh, indent=4)
    hist_mins_fig, hist_mins_axs = hist_minstops_per_seq(df)
    hist_mins_fig.savefig(args.out_hist_mins)
    dict_mins = mqc_hist_minstops_per_seq(df)
    with open(args.out_mqc_mins, "w") as fh:
        json.dump(dict_mins, fh, indent=4)
    if args.hmm:
        logger.info("Plotting histogram of HMM bitscores %s", args.out_hist_hmm)
        hist_hmm_fig, hist_hmm_axs = hist_hmm_scores(df, [llim, ulim])
        hist_hmm_fig.savefig(args.out_hist_hmm)
        dict_hmm_scores = mqc_hist_hmm_scores(df, [llim, ulim])
        with open(args.out_mqc_hmm, "w") as fh:
            json.dump(dict_hmm_scores, fh, indent=4)
        # Write screened sequences
        ok = list(df[df["ok"]]["seq_id"])
        logger.info("Writing ok sequences to file %s", args.out_screened)
        with open(args.out_screened, "w") as fh:
            SeqIO.write([nt[i] for i in ok], fh, "fasta")
        # logger.info("Writing bad sequences to file %s", args.out_bad)
        # with open(args.out_bad, "w") as fh:
        #     SeqIO.write([nt[i] for i in nt if i not in ok], fh, "fasta")
