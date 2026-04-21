from os import path
from numpy import isnan, isinf


class GCWriter:
    """
    This class is used for writing files like below:
    Sample  Chromosome    Start_position  End_position    value
    """
    def __init__(self, gc_file):
        self.gc_file = gc_file

    def write(self, gc_db):
        with open(self.gc_file, 'w') as fout:
            for _ in sorted(gc_db):
                chrn, sp, ep = _
                val = gc_db[_]
                fout.write("%s\t%d\t%d\t%d\n" % (chrn, sp, ep, val))


class BEDWriter:
    """
    This class is used for writing files like below:
    Sample  Chromosome    Start_position  End_position    value
    """
    def __init__(self, bed_file):
        self.bed_file = bed_file

    def write(self, bed_db):
        with open(self.bed_file, 'w') as fout:
            for smp in sorted(bed_db):
                for _ in sorted(bed_db[smp]):
                    chrn, sp, ep = _
                    val = bed_db[smp][_]
                    fout.write("%s\t%s\t%d\t%d\t%s\n" % (smp, chrn, sp, ep, str(val)))


class DepthWriter:
    """
    This class is used for writing sequence depth file which contain two columns:
    Sample  depth
    """
    def __init__(self, depth_file):
        self.depth_file = depth_file

    def write(self, depth_db):
        with open(self.depth_file, 'w') as fout:
            for smp in sorted(depth_db):
                fout.write("%s\t%f\n" % (smp, depth_db[smp]))


class GeneCNWriter:
    """
    This class is used for writing gc file like below:
    Gene    sample1_cn  sample2_cn  sample3_cn ...
    """
    def __init__(self, gene_cn_file):
        self.gene_cn_file = gene_cn_file

    def write(self, gene_cn, rounded=False):
        sample_list = []
        for _ in gene_cn:
            sample_list = sorted(gene_cn[_])
            break

        with open(self.gene_cn_file, 'w') as fout:
            fout.write("#Gene\t%s\n" % ('\t'.join(sample_list)))
            for gene in sorted(gene_cn):
                fout.write("%s" % gene)
                for smp in sample_list:
                    if smp not in gene_cn[gene]:
                        fout.write("\t%f" % float('nan'))
                    else:
                        if not rounded:
                            fout.write("\t%f" % gene_cn[gene][smp])
                        else:
                            if (not isnan(gene_cn[gene][smp])) and (not isinf(gene_cn[gene][smp])):
                                fout.write("\t%d" % int(gene_cn[gene][smp]))
                            else:
                                fout.write("\t%f" % gene_cn[gene][smp])
                fout.write("\n")


class TopRFDWriter:
    """
    This class is used for writing gc file like below:
    Gene    RFD p-value
    """
    def __init__(self, out_dir):
        self.out_dir = out_dir

    def write(self, rfd_db):
        gene_rfd_db = {}
        for grp in rfd_db:
            cnt = 0
            with open(path.join(self.out_dir, grp+".list"), 'w') as fout:
                for RFD, p, gn in sorted(rfd_db[grp], reverse=True):
                    if gn not in gene_rfd_db:
                        gene_rfd_db[gn] = {}
                    gene_rfd_db[gn][grp] = [RFD, p]
                    if cnt < len(rfd_db[grp])*.05 and p < .05:
                        fout.write("%s\t%f\t%f\n" % (gn, RFD, p))
                        cnt += 1

        with open(path.join(self.out_dir, "total.list"), 'w') as fout:
            fout.write("#Gene")
            for grp in sorted(rfd_db):
                fout.write("\t%s" % grp)
            fout.write("\n")

            for gn in sorted(gene_rfd_db):
                fout.write(gn)
                for grp in sorted(rfd_db):
                    if grp in gene_rfd_db[gn]:
                        fout.write("\t%f,% f" % (gene_rfd_db[gn][grp][0], gene_rfd_db[gn][grp][1]))
                    else:
                        fout.write("\tnan,nan")
                    fout.write("\n")
