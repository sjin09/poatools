from typing import Set, Dict, List, Tuple


class VCF:
    def __init__(self, line):
        arr = line.strip().split()
        self.chrom = arr[0]
        self.pos = int(arr[1])
        self.id = arr[2]
        self.ref = arr[3]
        self.alt_lst = arr[4].split(",")
        self.qual = arr[5]
        self.filter = True if arr[6] == "PASS" else False
        self.info = arr[7]
        self.format = arr[8]
        self.sample_format = arr[9]
        self.sample_format_lst = self.sample_format.split(":")
        self.sample_gt = self.sample_format_lst[0]
        if len(self.alt_lst) == 1: # bi-allelic
            self.is_biallelic = True
            self.alt = self.alt_lst[0]
            if len(self.ref) == 1 and len(self.alt) == 1: # snp
                self.is_snp = True
                self.is_indel = False
            else: # indel
                self.is_snp = False
                self.is_indel = True
        else:
            self.is_biallelic = False


def load_mut(
    vcf_file: str
) -> List[Tuple[str, int, str, str]]:
    mut_lst = []
    for i in open(vcf_file).readlines():
        if i.startswith("#"): continue
        v = VCF(i)
        if v.filter and v.is_snp and v.is_biallelic:
            mut_lst.append((v.chrom, v.pos, v.ref, v.alt))
    return mut_lst

