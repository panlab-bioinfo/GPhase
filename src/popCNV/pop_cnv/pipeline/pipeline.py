from os import makedirs, path, chdir
from gc import collect
from pop_cnv.io import message, loader, writer
from pop_cnv.worker import calculator
from pop_cnv.worker.dep_check import DepCheck


class Pipeline:
    def __init__(self, args):
        self.genome = path.abspath(args.genome)
        self.bam_path = path.abspath(args.bam)
        self.win_size = args.win_size
        self.gene_bed = path.abspath(args.list)
        self.workdir = path.abspath(args.workdir)
        self.group = path.abspath(args.group)
        self.wild_grp = args.wild
        self.threads = args.threads

    def run(self):
        msg = message.Message()
        step_list = ["00.Mosdepth",
                     "01.GC.stat",
                     "02.Sample.depth",
                     "03.rd",
                     "04.CN",
                     "05.genes.cn",
                     "06.genes.round.cn",
                     "07.RFD"]
        if not path.exists(self.workdir):
            makedirs(self.workdir)
        chdir(self.workdir)

        msg.info("Step00: Mosdepth")
        msg.info("Checking mosdepth...")
        mos_ava = DepCheck().check()
        if not mos_ava:
            msg.error("Mosdepth not found, exiting")
            exit(-1)
        mos_path = path.join(self.workdir, step_list[0])
        if not path.exists(mos_path):
            msg.info("Running mosdepth...")
            makedirs(mos_path)
            mos_runner = calculator.BamDepth()
            mos_runner.run(self.bam_path, self.win_size, mos_path, self.threads)
            msg.info("Done")
        else:
            msg.info("Directory %s found, skipping..." % mos_path)

        msg.info("Step01: GC Statistics")
        gc_file = path.join(self.workdir, step_list[1])
        if not path.exists(gc_file):
            msg.info("Calculating GC proportion...")
            gc_runner = calculator.GC()
            gc_runner.stat(self.genome, self.win_size, self.threads)
            gc_db = gc_runner.gc_db
            msg.info("Writing GC proportion...")
            gc_writer = writer.GCWriter(gc_file)
            gc_writer.write(gc_db)
            del gc_db
            collect()
            msg.info("Done")
        else:
            msg.info("File %s found, skipping..." % gc_file)

        msg.info("Step02: Depth of samples")
        depth_file = path.join(self.workdir, step_list[2])
        if not path.exists(depth_file):
            msg.info("Calculating depth")
            sd_runner = calculator.SeqDepth()
            sd_runner.calc(self.genome, self.bam_path, self.threads)
            sd_db = sd_runner.depth_db
            msg.info("Writing depth")
            sd_writer = writer.DepthWriter(depth_file)
            sd_writer.write(sd_db)
            msg.info("Done")
        else:
            msg.info("File %s found, loading..." % depth_file)
            sd_loader = loader.DepthLoader()
            sd_loader.load(depth_file)
            sd_db = sd_loader.depth_db
            msg.info("Loaded")

        msg.info("Step03: Read depth")
        rd_file = path.join(self.workdir, step_list[3])
        cn_file = path.join(self.workdir, step_list[4])
        rd_db = {}
        if not path.exists(rd_file):
            msg.info("Calculating read depth")
            rd_runner = calculator.Norm()
            rd_runner.norm(mos_path, sd_db, self.threads)
            rd_db = rd_runner.norm_db
            msg.info("Writing read depth")
            rd_writer = writer.BEDWriter(rd_file)
            rd_writer.write(rd_db)
            msg.info("Done")
        else:
            if not path.exists(cn_file):
                msg.info("File %s found, loading..." % rd_file)
                rd_loader = loader.BEDLoader()
                rd_loader.load(rd_file)
                rd_db = rd_loader.bed_db
                msg.info("Loaded")
            else:
                msg.info("File %s found, skipping..." % cn_file)

        msg.info("Step04: CN")
        gene_cn_file = path.join(self.workdir, step_list[5])
        cn_db = {}
        if not path.exists(cn_file):
            msg.info("Calculating CN")
            cn_runner = calculator.CN()
            cn_runner.convert(gc_file, rd_db, self.threads)
            del rd_db
            collect()
            cn_db = cn_runner.cn_db
            msg.info("Writing CN")
            cn_writer = writer.BEDWriter(cn_file)
            cn_writer.write(cn_db)
            msg.info("Done")
        else:
            if not path.exists(gene_cn_file):
                msg.info("File %s found, loading..." % cn_file)
                cn_loader = loader.BEDLoader()
                cn_loader.load(cn_file)
                cn_db = cn_loader.bed_db
                msg.info("Loaded")
            else:
                msg.info("File %s found, skipping..." % gene_cn_file)

        msg.info("Step05: Gene CN")
        round_cn_file = path.join(self.workdir, step_list[6])
        gene_cn_db = {}
        if not path.exists(gene_cn_file):
            msg.info("Calculating Gene CN")
            gene_loader = loader.GeneLoader()
            gene_loader.load(self.gene_bed)
            gene_bed_db = gene_loader.bed_db
            gene_cn_runner = calculator.GeneCN()
            gene_cn_runner.calc(cn_db, gene_bed_db)
            del cn_db
            collect()
            gene_cn_db = gene_cn_runner.gene_cn
            msg.info("Writing Gene CN")
            gene_cn_writer = writer.GeneCNWriter(gene_cn_file)
            gene_cn_writer.write(gene_cn_db)
            msg.info("Done")
        else:
            if not path.exists(round_cn_file):
                msg.info("File %s found, loading..." % gene_cn_file)
                gene_cn_loader = loader.GeneCNLoader()
                gene_cn_loader.load(gene_cn_file)
                gene_cn_db = gene_cn_loader.gene_cn
                msg.info("Loaded")
            else:
                msg.info("File %s found, skipping..." % round_cn_file)

        msg.info("Step06: Rounding CN")
        rfd_dir = path.join(self.workdir, step_list[7])
        round_cn_db = {}
        if not path.exists(round_cn_file):
            msg.info("Rounding Gene CN")
            round_cn_runner = calculator.RoundCN()
            round_cn_runner.round(gene_cn_db)
            del gene_cn_db
            collect()
            round_cn_db = round_cn_runner.round_cn
            msg.info("Writing Rounded Gene CN")
            round_cn_writer = writer.GeneCNWriter(round_cn_file)
            round_cn_writer.write(round_cn_db, rounded=True)
            msg.info("Done")
        else:
            if not path.exists(rfd_dir):
                msg.info("File %s found, loading..." % round_cn_file)
                round_cn_loader = loader.GeneCNLoader()
                round_cn_loader.load(round_cn_file)
                round_cn_db = round_cn_loader.gene_cn
                msg.info("Loaded")
            else:
                msg.info("Directory %s found, skipping..." % rfd_dir)

        msg.info("Step07: RFD and F-test")
        if not path.exists(rfd_dir):
            msg.info("Analyzing RFD")
            makedirs(rfd_dir)
            grp_loader = loader.GRPLoader()
            grp_loader.load(self.group)
            grp_db = grp_loader.grp_db
            rfd_runner = calculator.RFD()
            rfd_runner.calc(round_cn_db, grp_db, self.wild_grp)
            del round_cn_db
            collect()
            rfd_db = rfd_runner.rfd_db
            msg.info("Writing RFD")
            rfd_writer = writer.TopRFDWriter(rfd_dir)
            rfd_writer.write(rfd_db)
            msg.info("Done")
        else:
            msg.info("Directory %s found, completing" % rfd_dir)

        msg.info("All Completed")
