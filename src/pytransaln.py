#!/usr/bin/env python3

from Bio import SeqIO
from Bio import Align
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO
from subprocess import run, DEVNULL
from collections import defaultdict
import argparse

VERBOSE = True

# Future enhancements
# * Individually defined reading frames
# * User-supplied input amino acid alignment
# * Identify sequences with wrong frame or frameshifts


def translate_3_frames(seqdict, codes):
    out = defaultdict(lambda: defaultdict(int))
    for i in seqdict:
        for frame in [0,1,2]:
            newid = ";".join([i, f"frame={str(frame)}", f"code={str(codes[i])}"])
            out[i][frame] = seqdict[i][frame:].translate(table=codes[i], id=newid, name=newid)
    return out


def translate_1_frame(seqdict, frames, codes):
    out = defaultdict(lambda: defaultdict(int))
    for i in seqdict:
        newid = ";".join([i, f"frame={str(frames[i])}", f"code={str(codes[i])}"])
        out[i][frames[i]] = seqdict[i][frames[i]:].translate(table=codes[i], id=newid, name=newid)
    return out


def translate_minstops(seqdict, codes):
    threeframes = translate_3_frames(seqdict, codes)
    minstops = {}
    for i in threeframes:
        stopcounts = {frame : threeframes[i][frame].seq.count("*") for frame in [0,1,2]}
        minstops[i] = {frame : threeframes[i][frame] for frame in stopcounts if stopcounts[frame] == min(stopcounts.values())}
    return minstops


def guessframe(seqdict, codes):
    minstops = translate_minstops(seqdict, codes)
    # Assume that true reading frame has fewest stop codons
    ok = {i : minstops[i] for i in minstops if len(minstops[i]) == 1 }
    print(f"{str(len(ok))} of {str(len(minstops))} sequences have one frame with fewest stop codons")
    if len(ok) < len(minstops):
        print("Choosing reading frame for sequences with multiple minimal-stop frames by alignment scores")
        bestaln = {}
        aligner =  Align.PairwiseAligner()
        for i in minstops:
            if len(minstops[i])>1:
                alnscores = {frame : sum([aligner.score(minstops[i][frame], minstops[j][k]) for j in ok for k in ok[j]]) for frame in minstops[i]}
                bestframe = max(alnscores, key=lambda x: alnscores[x])
                bestaln[i] = {bestframe: minstops[i][bestframe]}
                if VERBOSE:
                    print(i)
                    print(alnscores)
        ok.update(bestaln)
    else:
        print("No ties to break")
    return ok


def trdict2seqlist(trdict):
    out = [trdict[i][frame] for i in trdict for frame in trdict[i]]
    aa2nt = {trdict[i][frame].id : i for i in trdict for frame in trdict[i]}
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
    for i in range(int(len(seq)/3)):
        yield seq[i*3 : i*3+3]
    if len(seq)%3 != 0: # trailing untranslated sequence
        yield seq[int(len(seq)/3)*3:]


def aa_aln_to_nt_aln(aa, nt, frame=0):
    """Align nucleotide sequence against aligned amino acid sequence
    
    Does not check for correctness of the translation. Unlike 
    Bio.codonalign, reading frame offset is taken into account.
    
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
    parser.add_argument("--input", help="Input unaligned nucleotide sequences, Fasta format")
    parser.add_argument("--guessframe", default=False, action="store_true", help="Guess best reading frame by minimizing stop codons")
    parser.add_argument("--frame", default=0, type=int, help="Reading frame offset to apply to all sequences, must be 0, 1, or 2")
    parser.add_argument("--code", default=5, type=int, help="Genetic code to use for all sequences, NCBI translation table number (except stopless codes 27, 28, 31)")
    parser.add_argument("--aligner", default="mafft", help="Alignment program to use (only MAFFT implemented at the moment)")
    parser.add_argument("--out_aa", default="test.aa.fasta", help="Path to write translated AA sequences")
    parser.add_argument("--out_aln_aa", default="test.aln.aa.fasta", help="Path to write aligned amino acid sequences")
    parser.add_argument("--out_aln_nt", default="test.aln.nt.fasta", help="Path to write aligned nucleotide sequences")
    parser.add_argument("--threads", default=1, type=int, help="Number of threads to pass to alignment program")
    args = parser.parse_args()

    if args.frame not in [0,1,2]:
        raise ValueError("Frame must be 0, 1, or 2")
        
    if args.code in [27, 28, 31]:
        raise ValueError("Please choose a non-ambiguous genetic code")

    nt = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))

    if guessframe:
        seq2code = { i : args.code for i in nt }
        tr = guessframe(seqdict=nt, codes=seq2code)
        seq2frame = { i : list(tr[i].keys())[0] for i in tr }
    else:
        # Single reading frame for all sequences
        seq2frame = { i : args.frame for i in nt }
        seq2code = { i : args.code for i in nt }
        tr = translate_1_frame(seqdict=nt, frames=seq2frame, codes=seq2code)

    aa, aa2nt = trdict2seqlist(tr)
    
    with open(args.out_aa, "w") as fh:
        SeqIO.write(aa, fh, "fasta")

    # read aa alignment
    print("Aligning with MAFFT")
    mafft_job = run(["mafft", "--thread", str(args.threads), args.out_aa], capture_output=True)
    print(mafft_job.stderr.decode())
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

if __name__ == "__main__":
    main()
