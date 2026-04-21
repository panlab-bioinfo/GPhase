from PySide6.QtCore import QObject, Signal, QSortFilterProxyModel
from os import path, makedirs
from pop_cnv.worker import calculator
from pop_cnv.io import loader, writer
from pop_cnv.ui.table_model import TableModel
from pop_cnv.worker import visualization
from gc import collect


class PopCNVWorker(QObject):
    progress = Signal(str)
    completed = Signal(int)

    def __init__(self):
        super(PopCNVWorker, self).__init__()
        self.__parent_form = None
        self.__genome_file = None
        self.__mosdepth_dir = None
        self.__seq_depth_file = None
        self.__smp_file = None
        self.__gene_list_file = None
        self.__win_size = None
        self.__wrkdir = None
        self.__wild_group = None
        self.__threads = None

    def set_param(self,
                  parent_form,
                  genome_file,
                  mosdepth_dir,
                  seq_depth_file,
                  smp_file,
                  gene_list_file,
                  win_size,
                  wrkdir,
                  wild_group,
                  threads):
        self.__parent_form = parent_form
        self.__genome_file = genome_file
        self.__mosdepth_dir = mosdepth_dir
        self.__seq_depth_file = seq_depth_file
        self.__smp_file = smp_file
        self.__gene_list_file = gene_list_file
        self.__win_size = win_size
        self.__wrkdir = wrkdir
        self.__wild_group = wild_group
        self.__threads = threads

    def run_popcnv(self, v):
        # Step 01
        self.progress.emit("GC statistics")

        gc_file = path.join(self.__wrkdir, "01.gc.txt")
        if not path.exists(gc_file):
            try:
                gc_runner = calculator.GC()
                gc_runner.stat(self.__genome_file, self.__win_size, self.__threads)
                gc_db = gc_runner.gc_db
                gc_writer = writer.GCWriter(gc_file)
                gc_writer.write(gc_db)
                del gc_db
                collect()
            except Exception as e:
                self.progress.emit(" Failed with: %s" % repr(e))
                return
        self.progress.emit("")

        # Step 02
        self.progress.emit("Loading sequence depth of samples")
        try:
            sd_loader = loader.DepthLoader()
            sd_loader.load(self.__seq_depth_file)
            sd_db = sd_loader.depth_db
        except Exception as e:
            self.progress.emit(" Failed with: %s" % repr(e))
            return
        self.progress.emit("")

        # Step 03
        rd_file = path.join(self.__wrkdir, "02.rd.txt")
        cn_file = path.join(self.__wrkdir, "03.cn.txt")
        rd_db = {}
        self.progress.emit("Calculating read depth")
        if not path.exists(rd_file):
            try:
                rd_runner = calculator.Norm()
                rd_runner.norm(self.__mosdepth_dir, sd_db, self.__threads)
                rd_db = rd_runner.norm_db
                rd_writer = writer.BEDWriter(rd_file)
                rd_writer.write(rd_db)
            except Exception as e:
                self.progress.emit(" Failed with: %s" % repr(e))
                return
        else:
            if not path.exists(cn_file):
                rd_loader = loader.BEDLoader()
                rd_loader.load(rd_file)
                rd_db = rd_loader.bed_db
        self.progress.emit("")

        # Step 04
        self.progress.emit("Calculating CN")
        gene_cn_file = path.join(self.__wrkdir, "04.gene_cn.txt")
        cn_db = {}
        if not path.exists(cn_file):
            try:
                cn_runner = calculator.CN()
                cn_runner.convert(gc_file, rd_db, self.__threads)
                del rd_db
                collect()
                cn_db = cn_runner.cn_db
                cn_writer = writer.BEDWriter(cn_file)
                cn_writer.write(cn_db)
            except Exception as e:
                self.progress.emit(" Failed with: %s" % repr(e))
                return
        else:
            if not path.exists(gene_cn_file):
                cn_loader = loader.BEDLoader()
                cn_loader.load(cn_file)
                cn_db = cn_loader.bed_db
        self.progress.emit("")

        # Step 05
        self.progress.emit("Calculating gene CN")
        round_cn_file = path.join(self.__wrkdir, "05.gene_cn_round.txt")
        gene_cn_db = {}
        if not path.exists(gene_cn_file):
            try:
                gene_loader = loader.GeneLoader()
                gene_loader.load(self.__gene_list_file)
                gene_bed_db = gene_loader.bed_db
                gene_cn_runner = calculator.GeneCN()
                gene_cn_runner.calc(cn_db, gene_bed_db)
                del cn_db
                collect()
                gene_cn_db = gene_cn_runner.gene_cn
                gene_cn_writer = writer.GeneCNWriter(gene_cn_file)
                gene_cn_writer.write(gene_cn_db)
            except Exception as e:
                self.progress.emit(" Failed with: %s" % repr(e))
                return
        else:
            if not path.exists(round_cn_file):
                gene_cn_loader = loader.GeneCNLoader()
                gene_cn_loader.load(gene_cn_file)
                gene_cn_db = gene_cn_loader.gene_cn
        self.progress.emit("")

        # Step 06
        self.progress.emit("Rounding CN")
        rfd_dir = path.join(self.__wrkdir, "06.RFD")
        round_cn_db = {}
        if not path.exists(round_cn_file):
            try:
                round_cn_runner = calculator.RoundCN()
                round_cn_runner.round(gene_cn_db)
                del gene_cn_db
                collect()
                round_cn_db = round_cn_runner.round_cn
                round_cn_writer = writer.GeneCNWriter(round_cn_file)
                round_cn_writer.write(round_cn_db, rounded=True)
            except Exception as e:
                self.progress.emit(" Failed with: %s" % repr(e))
                return
        else:
            if not path.exists(rfd_dir):
                round_cn_loader = loader.GeneCNLoader()
                round_cn_loader.load(round_cn_file)
                round_cn_db = round_cn_loader.gene_cn
        self.progress.emit("")

        # Step 07
        self.progress.emit("RFD and F-test")
        if not path.exists(rfd_dir):
            makedirs(rfd_dir)
            try:
                grp_loader = loader.GRPLoader()
                grp_loader.load(self.__smp_file)
                grp_db = grp_loader.grp_db
                rfd_runner = calculator.RFD()
                rfd_runner.calc(round_cn_db, grp_db, self.__wild_group)
                del round_cn_db
                collect()
                rfd_db = rfd_runner.rfd_db
                rfd_writer = writer.TopRFDWriter(rfd_dir)
                rfd_writer.write(rfd_db)
            except Exception as e:
                self.progress.emit(" Failed with: %s" % repr(e))
                return
        self.progress.emit("Done")
        self.completed.emit(1)


class DataPreviewWorker(QObject):
    progress = Signal(str)
    completed = Signal(int)

    def __init__(self):
        super(DataPreviewWorker, self).__init__()
        self.__parent_form = None
        self.__data_type = ""
        self.__data_file_path = ""
        self.__gene_list_file_path = ""
        self.__sample_name = ""
        self.__chr_name = ""
        self.__model = None
        self.table_data = []
        self.header_data = []

    def set_param(self,
                  parent_form,
                  data_type,
                  data_file_path,
                  gene_list_file_path,
                  sample_name,
                  chr_name):
        self.__parent_form = parent_form
        self.__data_type = data_type
        self.__data_file_path = data_file_path
        self.__gene_list_file_path = gene_list_file_path
        self.__sample_name = sample_name
        self.__chr_name = chr_name

    def run(self, v):
        self.__parent_form.ui.tablePreview.setModel(None)
        self.progress.emit("Loading %s" % self.__data_file_path)
        if path.isfile(self.__data_file_path):
            if "Gene" in self.__data_type:
                self.table_data = []
                self.header_data = ["Sample", "Gene", self.__data_type]
                gene_loader = loader.GeneLoader()
                gene_loader.load(self.__gene_list_file_path)
                data_loader = loader.GeneCNLoader()
                data_loader.load(self.__data_file_path)
                for gn in data_loader.gene_cn:
                    for smp in data_loader.gene_cn[gn]:
                        if self.__chr_name == "All" or self.__chr_name == gene_loader.bed_db[gn][0]:
                            if self.__sample_name == "All" or self.__sample_name == smp:
                                self.table_data.append([smp, gn, data_loader.gene_cn[gn][smp]])
                del data_loader.gene_cn
                del gene_loader.bed_db
                collect()
            else:
                self.table_data = []
                self.header_data = ["Sample", "Chromosome", "Start", "End", self.__data_type]
                bed_loader = loader.BEDLoader()
                bed_loader.load(self.__data_file_path)
                for smp in bed_loader.bed_db:
                    if self.__sample_name == "All" or self.__sample_name == smp:
                        for key in bed_loader.bed_db[smp]:
                            chrn, sp, ep = key
                            if self.__chr_name == "All" or self.__chr_name == chrn:
                                self.table_data.append([smp, chrn, sp, ep, bed_loader.bed_db[smp][key]])
                del bed_loader.bed_db
                collect()
        else:
            self.table_data = []
            self.header_data = ["Sample", "Gene", "RFD", "P-value"]
            gene_loader = loader.GeneLoader()
            gene_loader.load(self.__gene_list_file_path)

            if self.__sample_name == "All":
                full_fn = path.join(self.__data_file_path, "total.list")
                rfd_loader = loader.RFDLoader()
                rfd_loader.load(full_fn, file_type="all")
                for gn in rfd_loader.rfd_db:
                    if self.__chr_name == "All" or self.__chr_name == gene_loader.bed_db[gn][0]:
                        for smp in rfd_loader.rfd_db[gn]:
                            self.table_data.append([smp, gn,
                                                    rfd_loader.rfd_db[gn][smp][0], rfd_loader.rfd_db[gn][smp][1]])
                del rfd_loader.rfd_db
                collect()
            else:
                full_fn = path.join(self.__data_file_path, "%s.list" % self.__sample_name)
                rfd_loader = loader.RFDLoader()
                rfd_loader.load(full_fn, file_type="single")
                for gn in rfd_loader.rfd_db:
                    if self.__chr_name == "All" or self.__chr_name == gene_loader.bed_db[gn][0]:
                        self.table_data.append([self.__sample_name, gn,
                                                rfd_loader.rfd_db[gn][0], rfd_loader.rfd_db[gn][1]])
                del rfd_loader.rfd_db
                collect()
        self.__model = TableModel(self.table_data, self.header_data)
        proxymodel = QSortFilterProxyModel()
        proxymodel.setSourceModel(self.__model)

        self.__parent_form.ui.tablePreview.setModel(proxymodel)
        self.__parent_form.ui.tablePreview.resizeColumnsToContents()

        self.progress.emit("Done")
        self.completed.emit(1)


class PicDrawWorker(QObject):
    progress = Signal(str)
    completed = Signal(int)

    def __init__(self):
        super(PicDrawWorker, self).__init__()
        self.__parent_form = None
        self.__data_type = ""
        self.__data_file_path = ""
        self.__gene_list_file_path = ""
        self.__sample_name = ""
        self.__chr_name = ""
        self.__model = None
        self.__plot_data = {}

        self.pic = visualization.PlotContent()
        self.graph_scene = None

    def set_param(self,
                  parent_form,
                  data_type,
                  data_file_path,
                  gene_list_file_path,
                  sample_name,
                  chr_name):
        self.__parent_form = parent_form
        self.__data_type = data_type
        self.__data_file_path = data_file_path
        self.__gene_list_file_path = gene_list_file_path
        self.__sample_name = sample_name
        self.__chr_name = chr_name

    def run(self, v):
        stat = ""
        self.__plot_data = {}
        if path.isfile(self.__data_file_path):
            self.progress.emit("Loading %s" % self.__data_file_path)
            if "Gene" in self.__data_type:
                gene_loader = loader.GeneLoader()
                gene_loader.load(self.__gene_list_file_path)
                data_loader = loader.GeneCNLoader()
                data_loader.load(self.__data_file_path)
                for gn in data_loader.gene_cn:
                    for smp in data_loader.gene_cn[gn]:
                        chrn = gene_loader.bed_db[gn][0]
                        if self.__chr_name == "All" or self.__chr_name == chrn:
                            if self.__sample_name == "All" or self.__sample_name == smp:
                                if smp not in self.__plot_data:
                                    self.__plot_data[smp] = {}
                                if chrn not in self.__plot_data[smp]:
                                    self.__plot_data[smp][chrn] = {"X": [], "Y": []}
                                self.__plot_data[smp][chrn]["X"].append(gene_loader.bed_db[gn][1])
                                self.__plot_data[smp][chrn]["Y"].append(data_loader.gene_cn[gn][smp])
                del data_loader.gene_cn
                del gene_loader.bed_db
                collect()
            else:
                bed_loader = loader.BEDLoader()
                bed_loader.load(self.__data_file_path)
                for smp in bed_loader.bed_db:
                    if self.__sample_name == "All" or self.__sample_name == smp:
                        for key in bed_loader.bed_db[smp]:
                            chrn, sp, ep = key
                            if self.__chr_name == "All" or self.__chr_name == chrn:
                                if smp not in self.__plot_data:
                                    self.__plot_data[smp] = {}
                                if chrn not in self.__plot_data[smp]:
                                    self.__plot_data[smp][chrn] = {"X": [], "Y": []}
                                self.__plot_data[smp][chrn]["X"].append(sp)
                                self.__plot_data[smp][chrn]["Y"].append(bed_loader.bed_db[smp][key])
                del bed_loader.bed_db
                collect()
            stat = self.pic.gene_line_graph(self.__plot_data)
        else:
            gene_loader = loader.GeneLoader()
            gene_loader.load(self.__gene_list_file_path)

            full_fn = path.join(self.__data_file_path, "total.list")
            rfd_loader = loader.RFDLoader()
            rfd_loader.load(full_fn, file_type="all")
            for gn in rfd_loader.rfd_db:
                chrn = gene_loader.bed_db[gn][0]
                if self.__chr_name == "All" or self.__chr_name == chrn:
                    if chrn not in self.__plot_data:
                        self.__plot_data[chrn] = {"X": [], "Y": []}
                    self.__plot_data[chrn]["X"].append(gene_loader.bed_db[gn][1])
                    self.__plot_data[chrn]["Y"].append(rfd_loader.rfd_db[gn][self.__sample_name][1])
            del rfd_loader.rfd_db
            collect()
            stat = self.pic.gene_manhattan_graph(self.__plot_data)
        if stat == "Success":
            self.progress.emit("Done")
        else:
            self.progress.emit(stat)
        self.completed.emit(1)
