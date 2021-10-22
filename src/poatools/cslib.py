import os
import re
from typing import Set, Dict, List, Tuple

class BAM:
    def __init__(self, line):
        # target
        self.tname = line.reference_name
        self.tstart = int(line.reference_start)
        self.tend = int(line.reference_end)
        # query
        self.qname = line.query_name
        self.qstart = line.query_alignment_start
        self.qend = line.query_alignment_end
        self.zmw = self.qname.split("/")[1]
        self.qseq = line.query_sequence
        self.qlen = line.query_length
        self.bq_int_lst = line.query_qualities
        self.orientation = "-" if line.is_reverse else "+"
        hq_base_count = 0
        for i in self.bq_int_lst:
            if i == 93:
                hq_base_count +=1
        self.hq_base_fraction = hq_base_count/float(self.qlen)
        self.mapq = line.mapping_quality
        self.cs_tag = line.get_tag("cs") if line.has_tag("cs") else "."
        self.alignment_type = line.get_tag("tp") if line.has_tag("tp") else "."


def cs2lst(cs_tag):
    cslst = [cs for cs in re.split("(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)", cs_tag)]
    cslst = [cs.upper() for cs in cslst if cs != ""]
    return cslst


def cs2tuple(
    cs_str: str,
    qpos: int,
    qseq: str
) -> List[Tuple[int, str, str, int, int]]:

    cstuple_lst = []
    cs_lst = cs2lst(cs_str)
    for cs in cs_lst:
        m = cs[1:]
        mlen = len(m)
        qstart = qpos
        if cs.startswith("="):  # match # --cs=long
            cs = ":{}".format(mlen)
            t = (1, m, m, mlen, mlen)
        elif cs.startswith(":"):  # match # --cs=short
            mlen = int(m)
            qend = qpos + mlen
            m = qseq[qstart:qend]
            t = (1, m, m, mlen, mlen)
        elif cs.startswith("*"):  # snp # target and query
            mlen = 1
            ref, alt = list(m)
            t = (2, ref, alt, 1, 1)
        elif cs.startswith("+"):  # insertion # query
            ref = qseq[qpos - 1]
            alt = ref + m
            t = (3, ref, alt, 0, mlen)
        elif cs.startswith("-"):  # deletion # target
            alt = qseq[qpos - 1]
            ref = alt + m
            t = (4, ref, alt, mlen, 0)
            mlen = 0
        cstuple_lst.append(t)
        qpos += mlen
    return cstuple_lst


def blast_sequence_identity(cstuple_lst, qstart, qend, read_length):

    # soft-clipped bases
    umismatch = qstart  # soft-clipped bases
    dmismatch = read_length - qend
    soft_clip_length = umismatch + dmismatch

    # init
    match = 0
    length = soft_clip_length
    mismatch = soft_clip_length
    for cstuple in cstuple_lst:
        mstate, _ref, _alt, ref_len, alt_len = cstuple
        if mstate == 1:  # match
            length += ref_len
        elif mstate == 2:  # mismatch: snp
            length += 1
            mismatch += 1
        elif mstate == 3:  # mismatch: insertion
            length += alt_len
            mismatch += alt_len
        elif mstate == 4:  # mismatch: deletion
            length += ref_len
            mismatch += ref_len
    match = length - mismatch
    identity = match/float(length)
    return identity


def get_mut_state(tsub_lst, tindel_lst, tmbs_lst, tcomplex_lst):
    counts = [len(tsub_lst), len(tindel_lst), len(tmbs_lst), len(tcomplex_lst)]
    mut_state = [_count != 0 for _count in counts]
    return mut_state


def cs2mut(
    cstuple_lst: List[Tuple[int, str, str, int, int]],
    tname: str,
    tpos: int,
    zmw: str,
    qpos: int,
    min_bq: int,
    bq_int: List[int]
):
    # init: lst
    qbq_lst = []
    tsub_lst = []
    qsub_lst = []
    tmbs_lst = []
    tindel_lst = []
    tcomplex_lst = []
    # return: muts
    state = 0
    for cstuple in cstuple_lst:
        mstate, ref, alt, ref_len, alt_len, = cstuple
        if state == 0 and mstate == 1:  # init # match
            counter = 0
            state = 0
        elif state == 0 and mstate != 1:  # init # mismatch
            counter += 1
            ref_lst = [ref]
            alt_lst = [alt]
            if mstate == 2:  # snp
                state = 1
                tstart = tpos
                qstart = qpos
            elif mstate == 3 or mstate == 4:  # insertion # deletion
                state = 2
                tstart = tpos - 1
                qstart = qpos - 1
        elif state != 0 and mstate == 2:  # snp
            state = 1
            counter += 1
            ref_lst.append(ref)
            alt_lst.append(alt)
        elif state != 0 and mstate == 3:  # insertion
            state = 2
            counter += 1
            ref_lst.append(ref)
            alt_lst.append(alt)
        elif state != 0 and mstate == 4:  # deletion
            state = 2
            counter += 1
            ref_lst.append(ref)
            alt_lst.append(alt)
        elif (
            state != 0 and mstate == 1 and ref_len <= 10
        ):  # match # mnp: condition # snp, match, snp
            counter += 1
            state = state
            ref_lst.append(ref)
            alt_lst.append(ref)
        elif state != 0 and mstate == 1 and ref_len > 11:  # match # return
            state = 4
        tpos += ref_len # mismatch
        qpos += alt_len

        # return
        if state == 4:
            counts = len(ref_lst)
            ref = "".join(ref_lst)
            alt = "".join(alt_lst)
            ref_len = len(ref)
            alt_len = len(alt)
            tmut = (tname, tstart + 1, ref, alt)
            qmut = (zmw, qstart + 1, alt, ref)
            if counts == 1 and ref_len == 1 and alt_len == 1:  # snp
                tsub_lst.append(tmut)
                qsub_lst.append(qmut)
                qbq_lst.append(bq_int[qstart])
            elif counts == 1 and ref_len < alt_len:  # insertion
                tindel_lst.append(tmut)
            elif counts == 1 and ref_len > alt_len:  # deletion
                tindel_lst.append(tmut)
            elif counts > 1 and ref_len == alt_len:  # mbs
                tmbs_lst.append(tmut)
            elif counts > 1 and ref_len != alt_len:  # complex
                tcomplex_lst.append(tmut)
            state = 0  # init

    count = 0
    for bq in qbq_lst:
        if bq >= min_bq:
            count += 1
    bq_state = True if count > 0 else False
    mut_state = get_mut_state(tsub_lst, tindel_lst, tmbs_lst, tcomplex_lst)
    return mut_state, bq_state, tsub_lst, qsub_lst, qbq_lst, tmbs_lst, tindel_lst, tcomplex_lst

