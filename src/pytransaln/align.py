#!/usr/bin/env python3

from pytransaln.frameshifts import report_frameshifts
from pytransaln.translate import (
    translate_1_frame,
    onebestframe,
    guessframe,
)
from pytransaln.stats import summarize_framestats_with_hmm

import logging

import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO
from subprocess import run
from collections import defaultdict

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
    logger.info("%d nucleotide sequences to align", len(nt))
    meanlen = sum([len(i) for i in nt.values()]) / len(nt)
    logger.info("Mean sequence length %d", meanlen)
    if meanlen < 50:
        logger.warning(
            "Mean sequence length under 50 nt, stop codon based heuristic will be unreliable"
        )

    too_many_stops = {}
    if args.how.startswith("e"):
        logger.info(
            "Guessing reading frame for each sequence by minimizing stop codons"
        )
        seq2code = {i: args.code for i in nt}
        tr, too_many_stops = guessframe(
            seqdict=nt, codes=seq2code, maxstops=args.maxstops, ignore_terminal_stop=args.ignore_terminal_stop
        )
        seq2frame = {i: list(tr[i].keys())[0] for i in tr}
    elif args.how.startswith("c"):
        logger.info(
            "Find one reading frame for all sequence that minimizes total stop codons"
        )
        seq2code = {i: args.code for i in nt}
        tr, too_many_stops = onebestframe(
            seqdict=nt, codes=seq2code, maxstops=args.maxstops, ignore_terminal_stop=args.ignore_terminal_stop
        )
        seq2frame = {i: list(tr[i].keys())[0] for i in tr}
    else:
        logger.info(
            "Applying same reading frame offset %d for all sequences",
            args.frame,
        )
        # Single reading frame for all sequences
        seq2frame = {i: args.frame for i in nt}
        seq2code = {i: args.code for i in nt}
        tr, too_many_stops = translate_1_frame(
            seqdict=nt, frames=seq2frame, codes=seq2code, maxstops=args.maxstops, ignore_terminal_stop=args.ignore_terminal_stop
        )

    if len(too_many_stops) > 0:
        logger.info(
            "%d sequences with > %d stop codons",
            len(too_many_stops),
            args.maxstops,
        )
        if len(too_many_stops) >= 0.5 * (len(nt)):
            logger.warning(
                "More than 50% of sequences have too many stop codons; check genetic code and sequence orientation?"
            )

    if args.hmm:
        logger.info(
            "Using HMM %s to screen translations of %d sequences with <= %d stop codons",
            args.hmm,
            len(tr),
            args.maxstops,
        )
        df, llim, ulim = summarize_framestats_with_hmm(
            tr, args.hmm, args.out_hmmsearch, 1.5
        )
        # TODO output HMM screening stats
        ok = [(i.seq_id, i.frame) for i in df[df["hmm_ok"]].itertuples()]
        add = defaultdict(lambda: defaultdict(str))
        ok_dict = defaultdict(lambda: defaultdict(str))
        for i in tr:
            for frame in tr[i]:
                if (i, int(frame)) in ok:
                    ok_dict[i][frame] = tr[i][frame]
                else:
                    add[i][frame] = tr[i][frame]
        logger.info(
            "%d sequences failed HMM screening despite having <= %d stop codons; will not be used in initial alignment",
            len(add),
            args.maxstops,
        )
        too_many_stops.update(add)  # TODO separate these bad sequences in another file
        tr = ok_dict

    logger.info("%d sequences for initial alignment", len(tr))

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
            logger.info("Sequence %s has %d likely frameshifts", i, len(frameshifts[i]))
        dfs = {i: pd.DataFrame(frameshifts[i]) for i in frameshifts}
        for i in dfs:
            dfs[i]["seq_id"] = i
        df = pd.concat(list(dfs.values()))
        # Rearrange column order and set integer columns
        df = df[["seq_id", "refstart", "refend", "origstart", "origend", "gaplen"]]
        for intcol in ["refstart", "refend", "origstart", "origend", "gaplen"]:
            df[intcol] = df[intcol].astype(int)
        df.to_csv(args.out_bad_fs_report, sep="\t", index=False)
