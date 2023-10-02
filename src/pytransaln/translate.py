#!/usr/bin/env python3

import logging

from collections import defaultdict
from Bio.Align import PairwiseAligner

logger = logging.getLogger(__name__)


def translate_striptrailing(seq, frame, table, id, name):
    """Translate SeqRecord nucleotide sequence without trailing extra bases

    .translate() method of Bio.SeqRecord will complain if sequence length is
    not a multiple of three, ominously saying that this may be an error in
    future versions of Biopython. This function translates a SeqRecord with a
    specified frame offset and strips trailing extra bases so that the length
    is a multiple of three.

    Parameters
    ----------
    seq : Bio.SeqRecord
        Nucleotide sequence to be translated. Assumes no gaps else GIGO
    frame : int
        Reading frame offset, must be 0, 1, or 2 (forward) or -1, -2, -3 (reverse complement).
    table : int
        NCBI translation table number
    id : str
    name : str
        ID and Name for new SeqRecord object of the translation
    """
    if frame < 0:
        len_prime = len(seq) - (0 - frame - 1)
        trailing = (len_prime % 3)
        until = len(seq) - trailing
        return seq.reverse_complement()[-frame-1:until].translate(table=table, id=id, name=name)
    else:
        len_prime = len(seq) - frame
        trailing = (len_prime % 3)
        until = len(seq) - trailing
        return seq[frame:until].translate(table=table, id=id, name=name)


def translate_3_frames(seqdict, codes):
    """Translate nucleotide sequences into three forward frames

    Parameters
    ----------
    seqdict : dict
        SeqRecord objects for nucleotide sequences keyed by id
    codes : dict
        Genetic code for each sequence keyed by id

    Returns
    -------
    dict
        Translated sequences keyed by nucleotide sequence id (primary key) and
        frame offset (secondary key)
    """
    out = defaultdict(lambda: defaultdict(int))
    for i in seqdict:
        for frame in [0, 1, 2]:
            newid = ";".join([i, f"frame={str(frame)}", f"code={str(codes[i])}"])
            out[i][frame] = translate_striptrailing(
                seqdict[i], frame=frame, table=codes[i], id=newid, name=newid
            )
    return out


def onebestframe(seqdict, codes, maxstops):
    """Find one reading frame that minimizes stop codons for all sequences

    Assumes that all sequences have same frame, e.g. PCR amplicons with
    conserved primers.
    """
    trseq = translate_3_frames(seqdict, codes)
    sumstops = {
        frame: sum([trseq[i][frame].seq.count("*") for i in trseq])
        for frame in [0, 1, 2]
    }
    for frame in sumstops:
        logging.info("Frame %d has total %d stop codons", frame, sumstops[frame])
    bestframe = min(sumstops, key=lambda x: sumstops[x])
    ok = {
        i: {bestframe: trseq[i][bestframe]}
        for i in trseq
        if trseq[i][bestframe].seq.count("*") <= maxstops
    }
    too_many_stops = {
        i: {bestframe: trseq[i][bestframe]}
        for i in trseq
        if trseq[i][bestframe].seq.count("*") > maxstops
    }
    return ok, too_many_stops


def translate_1_frame(seqdict, frames, codes, maxstops):
    """Translate nucleotide sequences into a specified forward reading frame

    Parameters
    ----------
    seqdict : dict
        SeqRecord objects for nucleotide sequences keyed by id
    frames : dict
        Frame offset (0, 1, or 2) for each sequenced keyed by id
    codes : dict
        Genetic code for each sequence keyed by id

    Returns
    -------
    dict
        Translated sequences keyed by nucleotide sequence id (primary key) and
        frame offset (secondary key)
    """
    out = defaultdict(lambda: defaultdict(int))
    too_many_stops = defaultdict(lambda: defaultdict(int))
    for i in seqdict:
        newid = ";".join([i, f"frame={str(frames[i])}", f"code={str(codes[i])}"])
        trseq = translate_striptrailing(
            seqdict[i], frame=frames[i], table=codes[i], id=newid, name=newid
        )
        if trseq.seq.count("*") > maxstops:
            logger.info("%d stop codons in sequence %s", trseq.seq.count("*"), i)
            too_many_stops[i][frames[i]] = trseq
        else:
            out[i][frames[i]] = trseq
    return out, too_many_stops


def translate_minstops(seqdict, codes, maxstops):
    """Translate in all forward frames and report translation with fewest stops

    Parameters
    ----------
    seqdict : dict
        SeqRecord objects for nucleotide sequences keyed by id
    frames : dict
        Frame offset (0, 1, or 2) for each sequenced keyed by id
    maxstops : int
        Maximum number of stops to allow per sequence

    Returns
    -------
    (dict, dict)
        Tuple of two dicts. The first represents translated sequences keyed by
        nucleotide sequence id (primary key) and frame offset (secondary key),
        keeping frames with the fewest stop codons only (may be more than one
        with the same number), and with <= the max number of stop codons.  The
        second as above but containing sequences that have too many stop
        codons.
    """
    threeframes = translate_3_frames(seqdict, codes)
    minstops = {}
    too_many_stops = {}
    for i in threeframes:
        stopcounts = {
            frame: threeframes[i][frame].seq.count("*") for frame in [0, 1, 2]
        }
        if min(stopcounts.values()) <= maxstops:
            minstops[i] = {
                frame: threeframes[i][frame]
                for frame in stopcounts
                if stopcounts[frame] == min(stopcounts.values())
            }
        else:
            too_many_stops[i] = {
                frame: threeframes[i][frame]
                for frame in stopcounts
                if stopcounts[frame] == min(stopcounts.values())
            }
            logger.info(">= %d stop codons in sequence %s", min(stopcounts.values()), i)
    return minstops, too_many_stops


def guessframe(seqdict, codes, maxstops):
    """Translate and automatically find best reading frame offset

    For each nucleotide sequence, find reading frame that minimizes number of
    stop codons and where number of stop codons does not exceed maximum. If
    there is more than one frame with the same number of stop codons, then
    pairwise align each frame's translation to the translated "good" sequences,
    and pick the frame that maximizes alignment score.

    Sequences with too many stop codons are not included.

    Parameters
    ----------------------
    Same as translate_minstops
    """
    minstops, too_many_stops = translate_minstops(seqdict, codes, maxstops)
    # Assume that true reading frame has fewest stop codons
    ok = {i: minstops[i] for i in minstops if len(minstops[i]) == 1}
    logger.info(
        f"{str(len(ok))} of {str(len(minstops))} sequences have one frame with fewest stop codons"
    )
    if len(ok) < len(minstops):
        logger.info(
            "Choosing reading frame for sequences with multiple minimal-stop frames by alignment scores"
        )
        bestaln = {}
        aligner = PairwiseAligner()
        for i in minstops:
            if len(minstops[i]) > 1:
                alnscores = {
                    frame: sum(
                        [
                            aligner.score(minstops[i][frame], minstops[j][k])
                            for j in ok
                            for k in ok[j]
                        ]
                    )
                    for frame in minstops[i]
                }
                bestframe = max(alnscores, key=lambda x: alnscores[x])
                bestaln[i] = {bestframe: minstops[i][bestframe]}
                logger.info(i)
                logger.info(alnscores)
        ok.update(bestaln)
    else:
        logger.info("No ties to break")
    return ok, too_many_stops
