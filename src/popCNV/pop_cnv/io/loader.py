class GCLoader:
    """
    This class is used for loading gc file like below:
    Chromosome    Start_position  End_position    value
    """
    def __init__(self):
        self.bed_db = {}

    def load(self, gc_file):
        with open(gc_file, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                chrn = data[0]
                sp = int(data[1])
                ep = int(data[2])
                self.bed_db[tuple([chrn, sp, ep])] = int(data[3])


class BEDLoader:
    """
    This class is used for loading files like below:
    Sample  Chromosome    Start_position  End_position    value
    """
    def __init__(self):
        self.bed_db = {}

    def load(self, bed_file):
        with open(bed_file, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                smp = data[0]
                chrn = data[1]
                sp = int(data[2])
                ep = int(data[3])
                if smp not in self.bed_db:
                    self.bed_db[smp] = {}
                self.bed_db[smp][tuple([chrn, sp, ep])] = float(data[4])


class DepthLoader:
    """
    This class is used for loading sequence depth file which contain two columns:
    Sample  depth
    """
    def __init__(self):
        self.depth_db = {}

    def load(self, depth_file):
        with open(depth_file, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                self.depth_db[data[0]] = float(data[1])


class FastaLoader:
    def __init__(self):
        pass

    @staticmethod
    def load(fasta_file):
        fa_db = {}
        with open(fasta_file, 'r') as fin:
            sid = ''
            for line in fin:
                if line[0] == '>':
                    sid = line.strip()[1:]
                    fa_db[sid] = []
                else:
                    fa_db[sid].append(line.strip())
        for sid in fa_db:
            fa_db[sid] = ''.join(fa_db[sid])
        return fa_db


class GeneLoader:
    """
    This class is used for loading gc file like below:
    Chromosome    Start_position  End_position    GeneID
    """
    def __init__(self):
        self.bed_db = {}

    def load(self, gene_bed):
        with open(gene_bed, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                chrn = data[0]
                sp = int(data[1])
                ep = int(data[2])
                gn = data[3]
                self.bed_db[gn] = [chrn, sp, ep]


class GeneCNLoader:
    """
    This class is used for loading gc file like below:
    Gene    sample1_cn  sample2_cn  sample3_cn ...
    """
    def __init__(self):
        self.gene_cn = {}

    def load(self, gene_cn_file):
        with open(gene_cn_file, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                if line[0] == '#':
                    sample_list = data[1:]
                else:
                    gn = data[0]
                    self.gene_cn[gn] = {}
                    for _ in range(len(sample_list)):
                        self.gene_cn[gn][sample_list[_]] = float(data[_+1])


class GRPLoader:
    """
    This class is used for loading gc file like below:
    Sample1 Group1
    Sample2 Group1
    Sample3 Group2
    ...
    """
    def __init__(self):
        self.grp_db = {}

    def load(self, grp_file):
        with open(grp_file, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                self.grp_db[data[0]] = data[1]


class RFDLoader:
    def __init__(self):
        self.rfd_db = {}

    def load(self, rfd_file, file_type="all"):
        if file_type == "all":
            with open(rfd_file, 'r') as fin:
                for line in fin:
                    data = line.strip().split('\t')
                    if line[0] == '#':
                        smp_list = data[1:]
                    else:
                        gn = data[0]
                        self.rfd_db[gn] = {}
                        for _ in range(len(smp_list)):
                            rfd, pval = data[_+1].split(',')
                            self.rfd_db[gn][smp_list[_]] = [float(rfd), float(pval)]
        else:
            with open(rfd_file, 'r') as fin:
                for line in fin:
                    data = line.strip().split()
                    gn = data[0]
                    self.rfd_db[gn] = [float(data[1]), float(data[2])]
