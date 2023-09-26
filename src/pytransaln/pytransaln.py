#!/usr/bin/env python3

from Bio import SeqIO
from Bio import Align
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO
from subprocess import run
from collections import defaultdict
import argparse
import logging

logging.basicConfig(
    format="%(asctime)s : %(levelname)s : %(message)s", level=logging.DEBUG
)
logger = logging.getLogger(__name__)

# Future enhancements
# * Individually defined reading frames
# * User-supplied input amino acid alignment
# * Identify sequences with wrong frame or frameshifts
# * Identify likely frameshift positions from MAFFT .map file
# * Add pre and post frame sequence to alignment
# * Guess genetic code
# * Translate 6 frames


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
            out[i][frame] = seqdict[i][frame:].translate(
                table=codes[i], id=newid, name=newid
            )
    return out


def onebestframe(seqdict, codes, maxstops):
    """Find one reading frame that minimizes stop codons for all sequences

    Assumes that all sequences have same frame, e.g. PCR amplicons with
    conserved primers.
    """
    trseq = translate_3_frames(seqdict, codes)
    sumstops = { frame : sum([trseq[i][frame].seq.count("*") for i in trseq]) for frame in [0,1,2] }
    for frame in sumstops:
        logging.debug("Frame %d has total %d stop codons", frame, sumstops[frame])
    bestframe = min(sumstops, key=lambda x: sumstops[x])
    ok = { i : { bestframe : trseq[i][bestframe] } for i in trseq if trseq[i][bestframe].seq.count("*") <= maxstops }
    too_many_stops = { i : { bestframe : trseq[i][bestframe] } for i in trseq if trseq[i][bestframe].seq.count("*") > maxstops }
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
        trseq = seqdict[i][frames[i] :].translate(table=codes[i], id=newid, name=newid)
        if trseq.seq.count("*") > maxstops:
            logger.debug("%d stop codons in sequence %s", trseq.seq.count("*"), i)
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
            logger.debug(
                ">= %d stop codons in sequence %s", min(stopcounts.values()), i
            )
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
        aligner = Align.PairwiseAligner()
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
                logger.debug(i)
                logger.debug(alnscores)
        ok.update(bestaln)
    else:
        logger.info("No ties to break")
    return ok, too_many_stops


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", help="Input unaligned nucleotide sequences, Fasta format"
    )
    parser.add_argument(
        "--guessframe",
        default=False,
        action="store_true",
        help="Guess best reading frame for each sequence individually that minimizes stop codons; overrides --frame",
    )
    parser.add_argument(
        "--onebestframe",
        default=False,
        action="store_true",
        help="Find single reading frame for all sequences that minimizes stop codons; overrides --frame")
    parser.add_argument(
        "--maxstops",
        default=0,
        type=int,
        help="Max stop codons to allow in 'good' alignment; nt sequences over this threshold in all frames will be written to --out_bad",
    )
    parser.add_argument(
        "--frame",
        default=0,
        type=int,
        help="Reading frame offset to apply to all sequences, must be 0, 1, or 2",
    )
    parser.add_argument(
        "--code",
        default=5,
        type=int,
        help="Genetic code to use for all sequences, NCBI translation table number (except stopless codes 27, 28, 31)",
    )
    parser.add_argument(
        "--aligner",
        default="mafft",
        help="Alignment program to use (only MAFFT implemented at the moment)",
    )
    parser.add_argument(
        "--out_aa",
        default="test.aa.fasta",
        help="Path to write translated AA sequences",
    )
    parser.add_argument(
        "--out_bad",
        default="test.bad.nt.fasta",
        help="Path to write sequences with too many stop codons",
    )
    parser.add_argument(
        "--out_aln_aa",
        default="test.aln.aa.fasta",
        help="Path to write aligned amino acid sequences",
    )
    parser.add_argument(
        "--out_aln_nt",
        default="test.aln.nt.fasta",
        help="Path to write aligned nucleotide sequences",
    )
    parser.add_argument(
        "--out_aln_nt_aug",
        default="test.aln.nt.aug.fasta",
        help="Path to write aligned nucleotide sequences with likely frameshifted sequences added",
    )
    parser.add_argument(
        "--threads",
        default=1,
        type=int,
        help="Number of threads to pass to alignment program",
    )
    args = parser.parse_args()

    if args.frame not in [0, 1, 2]:
        raise ValueError("Frame must be 0, 1, or 2")

    if args.code in [27, 28, 31]:
        raise ValueError("Please choose a non-ambiguous genetic code")

    nt = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
    logger.info(f"{str(len(nt))} nucleotide sequences to align")

    too_many_stops = None
    if args.guessframe:
        logger.info(
            "Guessing reading frame for each sequence by minimizing stop codons"
        )
        seq2code = {i: args.code for i in nt}
        tr, too_many_stops = guessframe(
            seqdict=nt, codes=seq2code, maxstops=args.maxstops
        )
        seq2frame = {i: list(tr[i].keys())[0] for i in tr}
    elif args.onebestframe:
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
        if len(too_many_stops) >= 0.5*(len(nt)):
            logger.info("More than 50% of sequences have too many stop codons; check genetic code and sequence orientation?")
    logger.info(f"{str(len(tr))} sequences for initial alignment")

    aa, aa2nt = trdict2seqlist(tr)

    with open(args.out_aa, "w") as fh:
        SeqIO.write(aa, fh, "fasta")

    # read aa alignment
    logger.info("Aligning with MAFFT")
    cmd = ["mafft", "--thread", str(args.threads), args.out_aa]
    logger.debug("Command: %s", " ".join(cmd))
    mafft_job = run(cmd, capture_output=True)
    # logger.debug(mafft_job.stderr.decode())
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
        logger.debug("Command: %s", " ".join(cmd))
        mafft_add = run(
            cmd,
            capture_output=True,
        )
        # logger.debug(mafft_add.stderr.decode())
        with open(args.out_aln_nt_aug, "w") as fh:
            fh.write(mafft_add.stdout.decode())


if __name__ == "__main__":
    main()