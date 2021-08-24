import gzip
from Bio import SeqIO
from fastx.mutsig import purine, purine2pyrimidine, tri_lst


def fa_tricounts(infile: str, threshold: int, outfile):
    # init
    tri2counts = {_tri: 0 for _tri in tri_lst}
    fafile = (
        SeqIO.parse(infile, "fasta")
        if infile.endswith((".fa", ".fasta"))
        else SeqIO.parse(gzip.open(infile, "rt"), "fasta")
    )
    for fa in fafile:
        seq = str(fa.seq)
        slen = len(seq)
        if slen < threshold:
            continue
        seq_tri_lst = [seq[sdex : sdex + 3] for sdex, _base in enumerate(seq)]
        seq_tri_lst = [_ for _ in seq_tri_lst if len(_) == 3 and _.count("N") == 0]
        for tri in seq_tri_lst:
            if tri[1] in purine:
                tri = "".join([purine2pyrimidine[base] for base in tri[::-1]])
            tri2counts[tri] += 1
    for _tri in tri_lst:
        outfile.write("{}\t{}\n".format(_tri, tri2counts[_tri]))
    outfile.close()


def fq_tricounts(infile: str, outfile):
    # init
    tri2counts = {_tri: 0 for _tri in tri_lst}
    fqfile = open(infile) if infile.endswith((".fq", ".fastq")) else gzip.open(infile)

    for i, j in enumerate(fqfile):
        k = i % 4
        if k == 0:  # header
            continue
        elif k == 1:  # sequence
            seq = j.strip()
            seq = str(seq) if not isinstance(seq, bytes) else seq.decode("utf-8")
        elif k == 2:
            continue  # plus
        elif k == 3:  # quality
            seq_tri_lst = [seq[sdex : sdex + 3] for sdex, _base in enumerate(seq)]
            seq_tri_lst = [_ for _ in seq_tri_lst if len(_) == 3 and _.count("N") == 0]
            for tri in seq_tri_lst:
                if tri[1] in purine:
                    tri = "".join([purine2pyrimidine[base] for base in tri[::-1]])
                tri2counts[tri] += 1

    for _tri in tri_lst:
        outfile.write("{}\t{}\n".format(_tri, tri2counts[_tri]))


def tricounts(infile, threshold, outfile):
    if infile.endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz")):
        fa_tricounts(infile, threshold, outfile)
    elif infile.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz")):
        fq_tricounts(infile, outfile)
    else:
        print("tri_counts doesn't support the provided input")
