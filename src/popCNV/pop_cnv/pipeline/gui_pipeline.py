from os import path, listdir
from multiprocessing import cpu_count
from pop_cnv.ui import ui_popcnv_main_form
from pop_cnv.io import loader
from pop_cnv.ui.custom_control import ControlGraphicsScene
from pop_cnv.worker.gui_popcnv_worker import PopCNVWorker, DataPreviewWorker, PicDrawWorker
from PySide6.QtWidgets import QWidget, QFileDialog, QMessageBox
from PySide6.QtCore import QThread, Signal, QCoreApplication, QEventLoop
from traceback import format_exc
import gzip


class GPopCNV(QWidget):
    __run_signal = Signal(int)
    __preview_pic_signal = Signal(int)
    __preview_table_signal = Signal(int)

    def __init__(self):
        super(GPopCNV, self).__init__()
        self.ui = None

        # variables for popCNV
        self.__grp_db = {}
        self.__genome_file = ""
        self.__mosdepth_dir = ""
        self.__seq_depth_file = ""
        self.__smp_file = ""
        self.__gene_list_file = ""
        self.__win_size = 0
        self.__wrkdir = ""
        self.__threads = 0
        self.__chr_set = set()

        self.graph_scene = None

        # use QThread avoid main form freeze
        self.__popcnv_worker = PopCNVWorker()
        self.__popcnv_worker_thread = QThread()

        self.__draw_pic_worker = PicDrawWorker()
        self.__draw_pic_worker_thread = QThread()

        self.__data_preview_worker = DataPreviewWorker()
        self.__data_preview_worker_thread = QThread()

        self.__pic_export_filename = ""
        self.__table_export_filename = ""

        # data types
        self.__data_type = ["Read depth", "Copy number",
                            "Gene copy number", "Gene copy number (Rounded)", "RFD"]
        self.__data_file_db = {"Read depth": "02.rd.txt",
                               "Copy number": "03.cn.txt",
                               "Gene copy number": "04.gene_cn.txt",
                               "Gene copy number (Rounded)": "05.gene_cn_round.txt",
                               "RFD": "06.RFD"}
        self.__init_ui()

    def __init_ui(self):
        self.ui = ui_popcnv_main_form.Ui_frmMain()
        self.ui.setupUi(self)

        # Disable controls which cannot be used before run popCNV
        self.ui.grpSettings.setEnabled(False)
        self.ui.btn_run.setEnabled(False)
        self.ui.grpOutput.setEnabled(False)

        # Set button connects
        self.ui.btn_load_genome.clicked.connect(self.__get_genome_path)
        self.ui.btn_load_seqdepth.clicked.connect(self.__get_seqdepth_path)
        self.ui.btn_load_mosdepth.clicked.connect(self.__get_mosdepth_path)
        self.ui.btn_load_smp_list.clicked.connect(self.__get_smplist_path)
        self.ui.btn_load_gene_list.clicked.connect(self.__get_genelist_path)
        self.ui.btn_load_wrkdir.clicked.connect(self.__get_work_path)

        self.ui.btn_draw_pic.clicked.connect(self.__draw_pic)
        self.ui.btn_export_pic.clicked.connect(self.__export_pic)
        self.ui.btn_preview_table.clicked.connect(self.__preview_data)
        self.ui.btn_export_table.clicked.connect(self.__export_data)

        # Set combobox connects
        self.ui.cbox_table_data_type.currentTextChanged.connect(self.__modify_table_smp_cbox)
        self.ui.cbox_pic_data_type.currentTextChanged.connect(self.__modify_pic_smp_cbox)

        self.ui.btn_check.clicked.connect(self.__check_data)
        self.ui.btn_run.clicked.connect(self.__run_popcnv)

        # Enable sort tableview by click column header
        self.ui.tablePreview.setSortingEnabled(True)

    def __get_genome_path(self):
        genome_file_path = QFileDialog.getOpenFileName(self, "Load Genome")[0]
        if genome_file_path:
            self.ui.text_genome_file_path.setText(genome_file_path)

    def __get_mosdepth_path(self):
        mosdepth_path = QFileDialog.getExistingDirectory(self, "Load Mosdepth")
        if mosdepth_path:
            self.ui.text_mosdepth_folder_path.setText(mosdepth_path)

    def __get_seqdepth_path(self):
        seqdepth_path = QFileDialog.getOpenFileName(self, "Load Seq depth")[0]
        if seqdepth_path:
            self.ui.text_seq_depth_file_path.setText(seqdepth_path)

    def __get_smplist_path(self):
        smplist_path = QFileDialog.getOpenFileName(self, "Load sample list")[0]
        if smplist_path:
            self.ui.text_smp_list_file_path.setText(smplist_path)

    def __get_genelist_path(self):
        genelist_path = QFileDialog.getOpenFileName(self, "Load gene list")[0]
        if genelist_path:
            self.ui.text_gene_list_file_path.setText(genelist_path)

    def __get_work_path(self):
        work_path = QFileDialog.getExistingDirectory(self, "Select work directory")
        if work_path:
            self.ui.text_wrkdir.setText(work_path)

    def __modify_table_smp_cbox(self):
        self.ui.cbox_table_smp.clear()
        if self.ui.cbox_table_data_type.currentText() == "RFD":
            grp_list = ["All"]
            grp_set = set(self.__grp_db[_] for _ in self.__grp_db)
            for grp in sorted(grp_set):
                if grp != self.__wild_group:
                    grp_list.append(grp)
            self.ui.cbox_table_smp.addItems(grp_list)
        else:
            sample_list = ["All"]
            sample_list.extend(sorted(self.__grp_db))
            self.ui.cbox_table_smp.addItems(sample_list)

    def __modify_pic_smp_cbox(self):
        self.ui.cbox_pic_smp.clear()
        if self.ui.cbox_pic_data_type.currentText() == "RFD":
            grp_list = []
            grp_set = set(self.__grp_db[_] for _ in self.__grp_db)
            for grp in sorted(grp_set):
                if grp != self.__wild_group:
                    grp_list.append(grp)
            self.ui.cbox_pic_smp.addItems(grp_list)
        else:
            sample_list = ["All"]
            sample_list.extend(sorted(self.__grp_db))
            self.ui.cbox_pic_smp.addItems(sample_list)

    def __notify_with_title(self, info=""):
        if info:
            self.setWindowTitle("popCNV - %s" % info)
        else:
            self.setWindowTitle("popCNV")

    def __check_data(self):
        self.__genome_file = self.ui.text_genome_file_path.text()
        self.__mosdepth_dir = self.ui.text_mosdepth_folder_path.text()
        self.__seq_depth_file = self.ui.text_seq_depth_file_path.text()
        self.__smp_file = self.ui.text_smp_list_file_path.text()
        self.__gene_list_file = self.ui.text_gene_list_file_path.text()

        if not self.__genome_file or not path.isfile(self.__genome_file):
            QMessageBox.critical(self, "Loading error", "Genome file not found")
            return

        if not self.__mosdepth_dir or not path.isdir(self.__mosdepth_dir):
            QMessageBox.critical(self, "Loading error", "Mosdepth folder not found")
            return

        if not self.__seq_depth_file or not path.isfile(self.__seq_depth_file):
            QMessageBox.critical(self, "Loading error", "Seq depth file not found")
            return

        if not self.__smp_file or not path.isfile(self.__smp_file):
            QMessageBox.critical(self, "Loading error", "Sample list not found")
            return

        if not self.__gene_list_file or not path.isfile(self.__gene_list_file):
            QMessageBox.critical(self, "Loading error", "Gene list not found")
            return

        # Set threads list
        threads_list = [str(_) for _ in range(1, cpu_count() + 1)]
        self.ui.cbox_threads.clear()
        self.ui.cbox_threads.addItems(threads_list)
        self.ui.cbox_threads.setCurrentText(threads_list[-1])

        # Get chromosomes list
        self.__get_chr_ids()

        # Load sample list
        self.__notify_with_title("Loading samples")
        try:
            grp_loader = loader.GRPLoader()
            grp_loader.load(self.__smp_file)
            self.__grp_db = grp_loader.grp_db
            grp_set = set(self.__grp_db[_] for _ in self.__grp_db)
            self.ui.cbox_wild_group.clear()
            self.ui.cbox_wild_group.addItems(grp_set)
            self.__notify_with_title()
        except Exception as e:
            QMessageBox.critical(self, "Loading failed", format_exc())
            self.__notify_with_title(" Failed with: %s" % repr(e))
            return

        # Enable controls
        self.ui.grpSettings.setEnabled(True)
        self.ui.btn_run.setEnabled(True)

    def __get_win(self):
        for fn in listdir(self.__mosdepth_dir):
            if not fn.endswith(".regions.bed.gz"):
                continue
            mos_file = path.join(self.__mosdepth_dir, fn)
            with gzip.open(mos_file, 'rt') as fin:
                for line in fin:
                    data = line.strip().split()
                    sp = int(data[1])
                    ep = int(data[2])
                    return ep - sp

    def __get_chr_ids(self):
        with open(self.__gene_list_file) as fin:
            for line in fin:
                data = line.strip().split()
                self.__chr_set.add(data[0])

    def __popcnv_completed(self, v):
        self.__popcnv_worker_thread.quit()
        # Enable controls before popcnv finished
        self.ui.btn_run.setEnabled(True)
        self.ui.btn_check.setEnabled(True)
        self.ui.grpInput.setEnabled(True)
        self.ui.grpSettings.setEnabled(True)
        self.ui.grpOutput.setEnabled(True)

        # Set output options
        self.ui.cbox_pic_data_type.clear()
        self.ui.cbox_pic_data_type.addItems(self.__data_type)
        self.ui.cbox_table_data_type.clear()
        self.ui.cbox_table_data_type.addItems(self.__data_type)

        chr_list = ["All"]
        chr_list.extend(sorted(self.__chr_set))
        self.ui.cbox_pic_chr.clear()
        self.ui.cbox_pic_chr.addItems(chr_list)
        self.ui.cbox_table_chr.clear()
        self.ui.cbox_table_chr.addItems(chr_list)

    def __run_popcnv(self):
        # Check workdir available
        self.__wild_group = self.ui.cbox_wild_group.currentText()
        self.__wrkdir = self.ui.text_wrkdir.text()
        self.__threads = int(self.ui.cbox_threads.currentText())
        if not self.__wrkdir or not path.isdir(self.__wrkdir):
            QMessageBox.critical(self, "Running error", "Work directory not found")
            return

        # Disable controls before popcnv finished
        self.ui.btn_run.setEnabled(False)
        self.ui.btn_check.setEnabled(False)
        self.ui.grpInput.setEnabled(False)
        self.ui.grpSettings.setEnabled(False)

        self.__win_size = self.__get_win()
        self.__popcnv_worker.progress.connect(self.__notify_with_title)
        self.__popcnv_worker.completed.connect(self.__popcnv_completed)
        self.__popcnv_worker.set_param(self,
                                       self.__genome_file,
                                       self.__mosdepth_dir,
                                       self.__seq_depth_file,
                                       self.__smp_file,
                                       self.__gene_list_file,
                                       self.__win_size,
                                       self.__wrkdir,
                                       self.__wild_group,
                                       self.__threads)
        self.__run_signal.connect(self.__popcnv_worker.run_popcnv)
        self.__popcnv_worker.moveToThread(self.__popcnv_worker_thread)
        self.__run_signal.emit(1)
        self.__popcnv_worker_thread.start()

    def __data_preview_completed(self):
        self.__draw_pic_worker_thread.quit()
        # Enable controls
        self.ui.btn_run.setEnabled(True)
        self.ui.btn_check.setEnabled(True)
        self.ui.grpInput.setEnabled(True)
        self.ui.grpSettings.setEnabled(True)
        self.ui.btn_preview_table.setEnabled(True)
        self.ui.btn_export_table.setEnabled(True)

    def __preview_data(self):
        data_type = self.ui.cbox_table_data_type.currentText()
        sample_name = self.ui.cbox_table_smp.currentText()
        chr_name = self.ui.cbox_table_chr.currentText()
        data_file_path = path.join(self.__wrkdir, self.__data_file_db[data_type])

        self.__table_export_filename = "%s_%s_%s.csv" % (data_type, sample_name, chr_name)
        # Disable controls before popcnv finished
        self.ui.btn_run.setEnabled(False)
        self.ui.btn_check.setEnabled(False)
        self.ui.grpInput.setEnabled(False)
        self.ui.grpSettings.setEnabled(False)
        self.ui.btn_preview_table.setEnabled(False)
        self.ui.btn_export_table.setEnabled(False)

        self.__data_preview_worker.progress.connect(self.__notify_with_title)
        self.__data_preview_worker.completed.connect(self.__data_preview_completed)
        self.__data_preview_worker.set_param(self,
                                             data_type,
                                             data_file_path,
                                             self.__gene_list_file,
                                             sample_name,
                                             chr_name)
        self.__preview_table_signal.connect(self.__data_preview_worker.run)
        self.__data_preview_worker.moveToThread(self.__data_preview_worker_thread)
        self.__preview_table_signal.emit(1)
        self.__data_preview_worker_thread.start()

    def __export_data(self):
        export_file_path = QFileDialog.getSaveFileName(self, "Export Data",
                                                       self.__table_export_filename,
                                                       "csv(*.csv)")[0]
        if export_file_path and self.__data_preview_worker.table_data:
            try:
                with open(export_file_path, 'w') as fout:
                    fout.write("%s\n" % (','.join(self.__data_preview_worker.header_data)))
                    for info in self.__data_preview_worker.table_data:
                        fout.write("%s\n" % (','.join(map(str, info))))
                QMessageBox.information(self, "Export Data", "Exported Success")
            except Exception as e:
                QMessageBox.critical(self, "Export Data", "Failed to export")

    def __draw_pic_completed(self, v):
        self.__draw_pic_worker_thread.quit()
        # Draw picture
        if not self.graph_scene:
            self.graph_scene = ControlGraphicsScene()
            self.graph_scene.addWidget(self.__draw_pic_worker.pic.figure_content)
            self.ui.graphPreview.setScene(self.graph_scene)
            self.ui.graphPreview.show()
        else:
            self.__draw_pic_worker.pic.figure_content.draw()
        # Enable controls
        self.ui.btn_run.setEnabled(True)
        self.ui.btn_check.setEnabled(True)
        self.ui.grpInput.setEnabled(True)
        self.ui.grpSettings.setEnabled(True)
        self.ui.btn_draw_pic.setEnabled(True)
        self.ui.btn_export_pic.setEnabled(True)

    def __draw_pic(self):
        data_type = self.ui.cbox_pic_data_type.currentText()
        sample_name = self.ui.cbox_pic_smp.currentText()
        chr_name = self.ui.cbox_pic_chr.currentText()
        data_file_path = path.join(self.__wrkdir, self.__data_file_db[data_type])

        self.__pic_export_filename = "%s_%s_%s.pdf" % (data_type, sample_name, chr_name)
        # Disable controls before popcnv finished
        self.ui.btn_run.setEnabled(False)
        self.ui.btn_check.setEnabled(False)
        self.ui.grpInput.setEnabled(False)
        self.ui.grpSettings.setEnabled(False)
        self.ui.btn_draw_pic.setEnabled(False)
        self.ui.btn_export_pic.setEnabled(False)

        self.__draw_pic_worker.progress.connect(self.__notify_with_title)
        self.__draw_pic_worker.completed.connect(self.__draw_pic_completed)
        self.__draw_pic_worker.set_param(self,
                                         data_type,
                                         data_file_path,
                                         self.__gene_list_file,
                                         sample_name,
                                         chr_name)
        self.__preview_pic_signal.connect(self.__draw_pic_worker.run)
        self.__draw_pic_worker.moveToThread(self.__draw_pic_worker_thread)
        self.__preview_pic_signal.emit(1)
        self.__draw_pic_worker_thread.start()

    def __export_pic(self):
        export_file_path = QFileDialog.getSaveFileName(self, "Export Picture",
                                                       self.__pic_export_filename,
                                                       "pdf(*.pdf)")[0]
        if export_file_path and self.__draw_pic_worker.pic.figure_content:
            try:
                self.__draw_pic_worker.pic.figure_content.plt.savefig(export_file_path, bbox_inches='tight')
                QMessageBox.information(self, "Export Picture", "Exported Success")
            except Exception as e:
                QMessageBox.critical(self, "Export Picture", "Failed to export")
