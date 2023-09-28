from collections import defaultdict


def read_mafft_add_map(infile):
    """Read correspondence table from MAFFT --add alignment

    Table reported by MAFFT when --mapout option used. Coordinates are
    1-based inclusive.

    Returns
    -------
    dict
        keyed by sequence id of added sequence; secondary keys 'letter',
        'origpos', 'refpos' to list objects comprising the three columns of the
        MAFFT mapout table. Original sequence positions not represented in the
        reference alignment ("-" in ref alignment position column) are encoded
        here as integer -1.
    """
    out = defaultdict(lambda: defaultdict(list))
    currid = None
    with open(infile, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                next
            elif line.startswith(">"):
                id = line[1:].rstrip()
                currid = id
            else:
                [letter, origpos, refpos] = line.rstrip().split(", ")
                out[currid]["letter"].append(letter)
                out[currid]["origpos"].append(int(origpos))
                if refpos == "-":
                    out[currid]["refpos"].append(-1)
                else:
                    out[currid]["refpos"].append(int(refpos))
    return out


def contiguous_runs(seq, report="values"):
    """Report start/end value pairs of contiguous runs in an integer sequence

    Example: in the sequence 1,2,3,6,7,8,10,11,12, there are three runs:
    (1,3),(6,8),(10,12)

    Option report='values' returns the values from the sequence, while
    report='index' returns the indices instead
    """
    breakpoints = (
        [0] + [i for i in range(len(seq)) if seq[i] - seq[i - 1] > 1] + [len(seq)]
    )
    if report == "values":
        out = [
            (seq[breakpoints[i - 1]], seq[breakpoints[i] - 1])
            for i in range(1, len(breakpoints))
        ]
    else:
        out = [
            (breakpoints[i - 1], breakpoints[i] - 1) for i in range(1, len(breakpoints))
        ]
    return out


def report_inserts(refpos, origpos):
    """Report inserts relative to original sequence

    From correspondence of original sequence positions and reference alignment
    positions.  Coordinates are 1-based inclusive. Skip first and last three
    bases because they are prone to misalignment if there are long alignment
    gaps at the flanks.

    Returns
    -------
    list of dicts with coordinates and length of each gap
    """
    # from sequence of aligned positions, find gaps that are not multiples of three
    out = []
    for i in range(3, len(refpos) - 3):
        if refpos[i - 1] > 0 and refpos[i] - refpos[i - 1] > 1:
            # print("Gap of length " + str(refpos[i] - refpos[i-1] - 1) + " between original positions " + str(origpos[i]) + " " + str(origpos[i-1]))
            out.append(
                {
                    "refstart": refpos[i - 1],
                    "refend": refpos[i],
                    "origstart": origpos[i - 1],
                    "origend": origpos[i],
                    "gaplen": refpos[i] - refpos[i - 1] - 1,
                }
            )
    return out


def report_deletions(refpos, origpos):
    """Report deletions relative to original sequence

    Sequence coordinate in original not represented in reference alignment.
    Coordinates reported are 1-based inclusive.

    Returns
    -------
    list of dicts with coordinates and length of each gap
    """
    # TODO : ignore first run of gaps if alignment starts with gaps
    out = []
    gappos = [origpos[i] for i in range(len(refpos)) if refpos[i] == -1]
    if len(gappos) > 0:
        for a, b in contiguous_runs(gappos, report="values"):
            gaplen = b - a + 1
            # print("Gap of length " + str(gaplen) + " comprising original positions " + str(a) + " " + str(b) + " inclusive")
            out.append(
                {
                    "refstart": "-1",
                    "refend": "-1",
                    "origstart": a,
                    "origend": b,
                    "gaplen": gaplen,
                }
            )
    return out


def report_frameshifts(infile):
    """From MAFFT mapout table, report gaps that are potential frameshifts"""
    mappos = read_mafft_add_map(infile)
    out = {}
    for seq in mappos:
        inss = report_inserts(
            refpos=mappos[seq]["refpos"], origpos=mappos[seq]["origpos"]
        )
        dels = report_deletions(
            refpos=mappos[seq]["refpos"], origpos=mappos[seq]["origpos"]
        )
        out[seq] = [i for i in inss if i["gaplen"] % 3 != 0] + [
            i for i in dels if i["gaplen"] % 3 != 0
        ]
    return out
