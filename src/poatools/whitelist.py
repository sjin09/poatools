import os
import pysam
import poatools.cslib
import poatools.vcflib
from typing import Set, Dict, List, Tuple


def get_whitelist(
    bam_file: str,
    vcf_file: str,
    min_bq: int,
    min_mapq: int,
    min_hq_base_fraction: float,
    min_sequence_identity: float,
    whitelist_file: str
) -> None:

    cell2zmw = {}
    o = open(whitelist_file, "w")
    mut_lst = poatools.vcflib.load_mut(vcf_file)
    alignments = pysam.AlignmentFile(bam_file, "rb")
    o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("CHROM", "POS", "REF", "ALT", "ALT_COUNT", "SMRTcell", "ZMW"))
    for (chrom, pos, ref, alt) in mut_lst:
        sub_lst = []
        zmw_lst = []
        cell_lst = []
        orientation_lst = []
        alt_count = 0
        for j in alignments.fetch(chrom, pos, pos+1):
            read = poatools.cslib.BAM(j)
            if read.mapq < min_mapq or read.hq_base_fraction < min_hq_base_fraction or read.alignment_type != "P": # select primary alignments
                continue

            cstuple_lst = poatools.cslib.cs2tuple(read.cs_tag, read.qstart, read.qseq)
            sequence_identity = poatools.cslib.blast_sequence_identity(cstuple_lst, read.qstart, read.qend, read.qlen)
            if sequence_identity < min_sequence_identity:  # read filter
                continue

            mut_state, bq_state, tsub_lst, qsub_lst, _, _, _, _ = poatools.cslib.cs2mut(
                cstuple_lst,
                read.tname,
                read.tstart,
                read.zmw,
                read.qstart,
                min_bq,
                read.bq_int_lst,
            )
            if not bq_state or not mut_state[0]: # reads without substitutions
                continue

            if (chrom, pos, ref, alt) in tsub_lst:
                alt_count += 1
                cell, zmw = read.qname.split("/")[0:2]
                zmw_lst.append(zmw)
                cell_lst.append(cell)
                orientation_lst.append(read.orientation)
                qseq_bq = "".join([chr(bq + 33) for bq in read.bq_int_lst])
                if not os.path.exists(cell):
                    os.mkdir(cell)
                if cell not in cell2zmw:
                    cell2zmw[cell] = []
                cell2zmw[cell].append(read.zmw)
                _, qpos, _, _ = qsub_lst[tsub_lst.index((chrom, pos, ref, alt))]
                qsub = "{}:{}_{}/{}".format(zmw, qpos, alt, ref)
                sub_lst.append(qsub)
                fq = open(os.path.join(cell, "{}.fastq".format(read.zmw)), "w")
                fq.write("@{}\n{}\n+\n{}\n".format(read.qname, read.qseq, qseq_bq))
                fq.close()

        if alt_count == 1:
            o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, pos, ref, alt, alt_count, ",".join(cell_lst), ",".join(zmw_lst), ",".join(orientation_lst), ",".join(sub_lst)))

    for cell in cell2zmw:
        wlist = open("{}.whitelist".format(cell), "w")
        wlist.write("{}".format("\n".join(list(cell2zmw[cell]))))
        wlist.close()
    alignments.close()
    o.close()