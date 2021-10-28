import re
import os
import sys
import shutil
import pyfastx
import numpy as np 
from typing import List

def chunkstring(string, chunklen):
    chunks = [
        string[i : i + chunklen] for i in range(0, len(string), chunklen)
    ]
    return chunks


def load_fa(
    fa_file: str,
) -> List[str]:
    if fa_file.endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz")):
        seq_lst = []
        index_lst = []
        filtered_seq_lst = []
        fa = pyfastx.Fasta(fa_file, build_index=True)
        for seq in fa:
            seq_lst.append(str(seq))
        seq_len_lst = [len(seq) for seq in seq_lst]
        seq_median_len = np.median(seq_len_lst)
        for i, (seq, seq_len) in enumerate(zip(seq_lst, seq_len_lst)):
            if seq_len > 2*seq_median_len or seq_len < 0.5*seq_median_len:
                continue
            else:
                index_lst.append(i)
                filtered_seq_lst.append(seq)
        return index_lst, filtered_seq_lst
    else:
        sys.exit("Did you provide a FASTA file?")


def load_fq(
    fq_file: str
) -> List[str]:

    ccs = []
    for i, j in enumerate(open(fq_file).readlines()):
        k = i % 4
        if k == 0:  # header
            continue
        elif k == 1:  # sequence
            seq = j.strip().upper()
            seq = seq if not isinstance(seq, bytes) else seq.decode("utf-8")
        elif k == 2:
            continue  # plus
        elif k == 3:  # quality
            seq_bq = j.strip()
            seq_bq = seq_bq if not isinstance(seq_bq, bytes) else seq_bq.decode("utf-8")
            ccs.append((seq, seq_bq))
    return ccs[0]


def get_subread_poa(subread_file):
    subread_poa_file = subread_file.replace(".fasta", ".subreads.poa")
    os.system("abpoa {} -r 1 -s > {}".format(subread_file, subread_poa_file))


def get_ccs_subread_poa(cell, zmw, ccs_seq, subread_seq_lst):

    ccs_subread_seq = [ccs_seq] + subread_seq_lst
    ccs_subread_fasta_file = "./{}/poa/{}.msa.fasta".format(cell, zmw)
    ccs_subread_poa_file = ccs_subread_fasta_file.replace(".fasta", ".tmp.poa") 

    o = open(ccs_subread_fasta_file, "w")
    for i, seq in enumerate(ccs_subread_seq):
        if i == 0:
            o.write(">{}/ccs\n".format(zmw))
            for chunk in chunkstring(seq, 50):
                o.write("{}\n".format(chunk))
        else:
            o.write(">{}/subreads/{}\n".format(zmw, i))
            for chunk in chunkstring(seq, 50):
                o.write("{}\n".format(chunk))
    o.close()

    os.system("abpoa {} -r 1 -s > {}".format(ccs_subread_fasta_file, ccs_subread_poa_file))
    # os.system("rm {}".format(ccs_subread_fasta_file))
    return ccs_subread_poa_file


def load_msa(msa_file, ccs_seq_bq):
    count = 0
    ccs_msa = None
    subread_msa_lst = []
    for line in open(msa_file).readlines():
        if line.startswith(">"): 
            continue
        if count == 0:
            count += 1
            ccs_msa = line.strip()
        else:
            count += 1
            subread_msa_lst.append(line.strip())
    # os.system("rm {}".format(msa_file))

    if ccs_msa is None:
        sys.exit("partial order alignment failed")

    index = 0
    ccs_msa_bq_lst = []
    for ccs_msa_base in ccs_msa:
        if ccs_msa_base == "-":
            ccs_msa_bq_lst.append("!")
        else:
            ccs_msa_bq_lst.append(ccs_seq_bq[index])
            index += 1
    ccs_msa_bq = "".join(ccs_msa_bq_lst)
    return ccs_msa, ccs_msa_bq, subread_msa_lst


def get_msa_pos2ccs_pos(ccs_msa):
    idx = 0
    msa2ccs = {}
    for jdx, ccs_msa_base in enumerate(ccs_msa):
        if ccs_msa_base == "-": 
            continue
        msa2ccs[jdx] = idx
        idx += 1
    return msa2ccs


def qsub_split(qsub):
    zmw, qpos, alt, ref = re.split(":|_|/", qsub)
    qpos = int(qpos) - 1
    return zmw, qpos, alt, ref


def get_subset_poa_pdf(cell, qsub, x, ccs_msa_chunk_lst, subread_msa_chunk_lst):

    zmw, qpos, alt, ref = qsub_split(qsub)
    subset_pdf_file = "./{}/poa/{}_{}_{}_{}.pdf".format(cell, zmw, qpos, alt, ref)
    subset_fasta_file = subset_pdf_file.replace(".pdf", ".fasta")

    fa = open(subset_fasta_file, "w")
    fa.write(">{}/ccs\n{}\n".format(zmw, "".join([_ for _ in ccs_msa_chunk_lst[x] if _ != "-"])))
    for q, subread_chunk_lst in enumerate(subread_msa_chunk_lst):
        subread_subset = "".join([_ for _ in subread_chunk_lst[x] if _ != "-"])
        if len(subread_subset) != 0:
            fa.write(">{}/subread/{}\n{}\n".format(zmw, q, subread_subset))
    fa.close()
    os.system("abpoa -g {} -r3 {} > /dev/null".format(subset_pdf_file, subset_fasta_file))
    os.system("rm {}.dot".format(subset_pdf_file)) 
    os.system("rm {}".format(subset_fasta_file))



def get_poa(cell, zmw, ccs_file, subread_file, qsub_lst, mlen, poa_file, pdf_state):

    if shutil.which("abpoa") is None:
        sys.exit("abPOA is not installed on your system")

    # load
    # get_subread_poa(subread_file) # subread poa
    qsub_lst = qsub_lst.split(",")
    ccs_seq, ccs_seq_bq = load_fq(ccs_file)
    subread_index_lst, subread_seq_lst = load_fa(subread_file)
    ccs_subread_poa_file = get_ccs_subread_poa(cell, zmw, ccs_seq, subread_seq_lst)
    ccs_msa, ccs_msa_bq, subread_msa_lst = load_msa(ccs_subread_poa_file, ccs_seq_bq)
    ccs_msa_chunk_lst = chunkstring(ccs_msa, mlen)
    ccs_msa_bq_chunk_lst = chunkstring(ccs_msa_bq, mlen)
    subread_msa_chunk_lst = [chunkstring(_, mlen) for _ in subread_msa_lst]
    chunk_len = len(ccs_msa_chunk_lst) - 1
    msa_pos2ccs_pos = get_msa_pos2ccs_pos(ccs_msa)

    # main
    counter = 0 
    poa = open(poa_file, "w")
    poa_dir = "/".join(poa_file.split("/")[:-1])
    for qsub in qsub_lst:
        _, qpos, alt, ref = qsub_split(qsub)
        poa_subset_file = os.path.join(poa_dir, "{}_{}_{}_{}.poa".format(zmw, qpos, ref, alt))
        poa_subset = open(poa_subset_file, "w")
        for x, ccs_msa_bq_chunk in enumerate(ccs_msa_bq_chunk_lst):
            if ccs_msa_bq_chunk.count("!") == len(ccs_msa_bq_chunk):
                continue
            for y, ccs_msa_bq in enumerate(ccs_msa_bq_chunk):
                if ccs_msa_bq != "!": 
                    break
            ccs_msa_start = (x * mlen) + y

            if x != chunk_len:
                ccs_msa_bq_chunk_reverse = ccs_msa_bq_chunk[::-1]
                for z, ccs_msa_bq in enumerate(ccs_msa_bq_chunk_reverse):
                    if ccs_msa_bq != "!": 
                        break
                ccs_msa_end = ((x * mlen) + mlen - 1) - z
            else:
                ccs_msa_bq_chunk_reverse = ccs_msa_bq_chunk[::-1]
                for z, ccs_msa_bq in enumerate(ccs_msa_bq_chunk_reverse):
                    if ccs_msa_bq != "!": 
                        break
                ccs_msa_end = ((x * mlen) + len(ccs_msa_bq_chunk_reverse) - 1) - z
            ccs_base_start = msa_pos2ccs_pos[ccs_msa_start]
            ccs_base_end = msa_pos2ccs_pos[ccs_msa_end]

            if qpos >= ccs_base_start and qpos <= ccs_base_end:
                ccs_msa_base_pos = (x * mlen)
                ccs_msa_chunk = ccs_msa_chunk_lst[x]
                for p, ccs_msa_base in enumerate(ccs_msa_chunk):
                    if not ccs_msa_base == "-":
                        ccs_base_pos =  msa_pos2ccs_pos[ccs_msa_base_pos]
                        if ccs_base_pos == qpos:
                            qsub_watermark_index = p
                            break
                    ccs_msa_base_pos += 1
                qsub_watermark = list(" " * mlen)
                qsub_watermark[qsub_watermark_index-1:qsub_watermark_index+2] = ["<","*",">"]
                qsub_watermark = "".join(qsub_watermark)

                msa_chunk_lst = []
                msa_chunk_lst.append("{:<15s}\t{}".format("qsub:{}_{}/{}".format(qpos, alt, ref), qsub_watermark))
                msa_chunk_lst.append("{:<15s}\t{}".format("bq:", ccs_msa_bq_chunk)) 
                msa_chunk_lst.append("{:<15s}\t{}".format("ccs:{}-{}".format(ccs_base_start, ccs_base_end), ccs_msa_chunk_lst[x]))
                for q, subread_chunk_lst in enumerate(subread_msa_chunk_lst):
                    msa_chunk_lst.append("{:<15s}\t{}".format("subread:{}".format(subread_index_lst[q]), subread_chunk_lst[x]))
                poa_subset.write("{}".format("\n".join(msa_chunk_lst + [""])))
                poa_subset.close()
                if counter == 0:
                    poa.write("{}\n".format("\n".join(msa_chunk_lst + [""])))
                if pdf_state:
                    get_subset_poa_pdf(cell, qsub, x, ccs_msa_chunk_lst, subread_msa_chunk_lst)
            else:
                msa_chunk_lst = []
                msa_chunk_lst.append("{:<15s}\t{}".format("bq:", ccs_msa_bq_chunk)) 
                msa_chunk_lst.append("{:<15s}\t{}".format("ccs:{}-{}".format(ccs_base_start, ccs_base_end), ccs_msa_chunk_lst[x]))
                for q, subread_chunk_lst in enumerate(subread_msa_chunk_lst):
                    msa_chunk_lst.append("{:<15s}\t{}".format("subread:{}".format(subread_index_lst[q]), subread_chunk_lst[x]))
                if counter == 0:
                    poa.write("{}\n".format("\n".join(msa_chunk_lst + [""])))
        counter += 1
        poa.close()

