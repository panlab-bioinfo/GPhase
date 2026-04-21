# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'popcnv_main_formFEiWUd.ui'
##
## Created by: Qt User Interface Compiler version 6.6.2
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
                            QMetaObject, QObject, QPoint, QRect,
                            QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QBrush, QColor, QConicalGradient, QCursor,
                           QFont, QFontDatabase, QGradient, QIcon,
                           QImage, QKeySequence, QLinearGradient, QPainter,
                           QPalette, QPixmap, QRadialGradient, QTransform)
from PySide6.QtWidgets import (QApplication, QComboBox, QGraphicsView, QGridLayout,
                               QGroupBox, QHBoxLayout, QHeaderView, QLabel,
                               QLineEdit, QPushButton, QSizePolicy, QSpacerItem,
                               QTabWidget, QTableView, QVBoxLayout, QWidget)
from pop_cnv.ui.custom_control import DragLineEdit, ControlGraphicsView


class Ui_frmMain(object):
    def setupUi(self, frmMain):
        if not frmMain.objectName():
            frmMain.setObjectName(u"frmMain")
        frmMain.resize(960, 1064)
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(frmMain.sizePolicy().hasHeightForWidth())
        frmMain.setSizePolicy(sizePolicy)
        self.gridLayout_7 = QGridLayout(frmMain)
        self.gridLayout_7.setObjectName(u"gridLayout_7")
        self.gridLayout_6 = QGridLayout()
        self.gridLayout_6.setObjectName(u"gridLayout_6")
        self.grpSettings = QGroupBox(frmMain)
        self.grpSettings.setObjectName(u"grpSettings")
        font = QFont()
        font.setFamilies([u"JetBrains Mono"])
        font.setPointSize(14)
        font.setBold(True)
        self.grpSettings.setFont(font)
        self.gridLayout_2 = QGridLayout(self.grpSettings)
        self.gridLayout_2.setObjectName(u"gridLayout_2")
        self.verticalLayout_2 = QVBoxLayout()
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.horizontalLayout_7 = QHBoxLayout()
        self.horizontalLayout_7.setObjectName(u"horizontalLayout_7")
        self.label_5 = QLabel(self.grpSettings)
        self.label_5.setObjectName(u"label_5")
        font1 = QFont()
        font1.setFamilies([u"JetBrains Mono"])
        font1.setPointSize(12)
        font1.setBold(True)
        self.label_5.setFont(font1)

        self.horizontalLayout_7.addWidget(self.label_5)

        self.cbox_wild_group = QComboBox(self.grpSettings)
        self.cbox_wild_group.setObjectName(u"cbox_wild_group")
        font2 = QFont()
        font2.setFamilies([u"JetBrains Mono"])
        font2.setPointSize(14)
        font2.setBold(False)
        self.cbox_wild_group.setFont(font2)

        self.horizontalLayout_7.addWidget(self.cbox_wild_group)

        self.horizontalSpacer = QSpacerItem(98, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_7.addItem(self.horizontalSpacer)

        self.label_6 = QLabel(self.grpSettings)
        self.label_6.setObjectName(u"label_6")
        self.label_6.setFont(font1)

        self.horizontalLayout_7.addWidget(self.label_6)

        self.cbox_threads = QComboBox(self.grpSettings)
        self.cbox_threads.setObjectName(u"cbox_threads")
        self.cbox_threads.setFont(font2)

        self.horizontalLayout_7.addWidget(self.cbox_threads)

        self.verticalLayout_2.addLayout(self.horizontalLayout_7)

        self.horizontalLayout_6 = QHBoxLayout()
        self.horizontalLayout_6.setObjectName(u"horizontalLayout_6")
        self.label_7 = QLabel(self.grpSettings)
        self.label_7.setObjectName(u"label_7")
        self.label_7.setFont(font1)

        self.horizontalLayout_6.addWidget(self.label_7)

        self.text_wrkdir = DragLineEdit(self.grpSettings)
        self.text_wrkdir.setObjectName(u"text_wrkdir")
        self.text_wrkdir.setFont(font2)

        self.horizontalLayout_6.addWidget(self.text_wrkdir)

        self.btn_load_wrkdir = QPushButton(self.grpSettings)
        self.btn_load_wrkdir.setObjectName(u"btn_load_wrkdir")
        self.btn_load_wrkdir.setFont(font2)

        self.horizontalLayout_6.addWidget(self.btn_load_wrkdir)

        self.verticalLayout_2.addLayout(self.horizontalLayout_6)

        self.gridLayout_2.addLayout(self.verticalLayout_2, 0, 0, 1, 1)

        self.gridLayout_6.addWidget(self.grpSettings, 2, 0, 1, 1)

        self.grpOutput = QGroupBox(frmMain)
        self.grpOutput.setObjectName(u"grpOutput")
        self.grpOutput.setFont(font)
        self.gridLayout_5 = QGridLayout(self.grpOutput)
        self.gridLayout_5.setObjectName(u"gridLayout_5")
        self.tabOutput = QTabWidget(self.grpOutput)
        self.tabOutput.setObjectName(u"tabOutput")
        self.tabOutput.setFont(font)
        self.tab_pic = QWidget()
        self.tab_pic.setObjectName(u"tab_pic")
        self.gridLayout_3 = QGridLayout(self.tab_pic)
        self.gridLayout_3.setObjectName(u"gridLayout_3")
        self.verticalLayout_3 = QVBoxLayout()
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.graphPreview = ControlGraphicsView(self.tab_pic)
        self.graphPreview.setObjectName(u"graphPreview")

        self.verticalLayout_3.addWidget(self.graphPreview)

        self.horizontalLayout_8 = QHBoxLayout()
        self.horizontalLayout_8.setObjectName(u"horizontalLayout_8")
        self.label_9 = QLabel(self.tab_pic)
        self.label_9.setObjectName(u"label_9")
        self.label_9.setFont(font1)

        self.horizontalLayout_8.addWidget(self.label_9)

        self.cbox_pic_data_type = QComboBox(self.tab_pic)
        self.cbox_pic_data_type.setObjectName(u"cbox_pic_data_type")
        sizePolicy1 = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.cbox_pic_data_type.sizePolicy().hasHeightForWidth())
        self.cbox_pic_data_type.setSizePolicy(sizePolicy1)
        self.cbox_pic_data_type.setFont(font2)

        self.horizontalLayout_8.addWidget(self.cbox_pic_data_type)

        self.label_10 = QLabel(self.tab_pic)
        self.label_10.setObjectName(u"label_10")
        self.label_10.setFont(font1)

        self.horizontalLayout_8.addWidget(self.label_10)

        self.cbox_pic_smp = QComboBox(self.tab_pic)
        self.cbox_pic_smp.setObjectName(u"cbox_pic_smp")
        sizePolicy1.setHeightForWidth(self.cbox_pic_smp.sizePolicy().hasHeightForWidth())
        self.cbox_pic_smp.setSizePolicy(sizePolicy1)
        self.cbox_pic_smp.setFont(font2)

        self.horizontalLayout_8.addWidget(self.cbox_pic_smp)

        self.label_11 = QLabel(self.tab_pic)
        self.label_11.setObjectName(u"label_11")
        self.label_11.setFont(font1)

        self.horizontalLayout_8.addWidget(self.label_11)

        self.cbox_pic_chr = QComboBox(self.tab_pic)
        self.cbox_pic_chr.setObjectName(u"cbox_pic_chr")
        sizePolicy1.setHeightForWidth(self.cbox_pic_chr.sizePolicy().hasHeightForWidth())
        self.cbox_pic_chr.setSizePolicy(sizePolicy1)
        self.cbox_pic_chr.setFont(font2)

        self.horizontalLayout_8.addWidget(self.cbox_pic_chr)

        self.horizontalSpacer_2 = QSpacerItem(118, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_8.addItem(self.horizontalSpacer_2)

        self.btn_draw_pic = QPushButton(self.tab_pic)
        self.btn_draw_pic.setObjectName(u"btn_draw_pic")
        self.btn_draw_pic.setFont(font2)

        self.horizontalLayout_8.addWidget(self.btn_draw_pic)

        self.btn_export_pic = QPushButton(self.tab_pic)
        self.btn_export_pic.setObjectName(u"btn_export_pic")
        self.btn_export_pic.setFont(font2)

        self.horizontalLayout_8.addWidget(self.btn_export_pic)

        self.verticalLayout_3.addLayout(self.horizontalLayout_8)

        self.gridLayout_3.addLayout(self.verticalLayout_3, 0, 0, 1, 1)

        self.tabOutput.addTab(self.tab_pic, "")
        self.tab_table = QWidget()
        self.tab_table.setObjectName(u"tab_table")
        self.gridLayout_4 = QGridLayout(self.tab_table)
        self.gridLayout_4.setObjectName(u"gridLayout_4")
        self.verticalLayout_4 = QVBoxLayout()
        self.verticalLayout_4.setObjectName(u"verticalLayout_4")
        self.tablePreview = QTableView(self.tab_table)
        self.tablePreview.setObjectName(u"tablePreview")

        self.verticalLayout_4.addWidget(self.tablePreview)

        self.horizontalLayout_9 = QHBoxLayout()
        self.horizontalLayout_9.setObjectName(u"horizontalLayout_9")
        self.label_12 = QLabel(self.tab_table)
        self.label_12.setObjectName(u"label_12")
        self.label_12.setFont(font1)

        self.horizontalLayout_9.addWidget(self.label_12)

        self.cbox_table_data_type = QComboBox(self.tab_table)
        self.cbox_table_data_type.setObjectName(u"cbox_table_data_type")
        sizePolicy1.setHeightForWidth(self.cbox_table_data_type.sizePolicy().hasHeightForWidth())
        self.cbox_table_data_type.setSizePolicy(sizePolicy1)
        self.cbox_table_data_type.setFont(font2)

        self.horizontalLayout_9.addWidget(self.cbox_table_data_type)

        self.label_14 = QLabel(self.tab_table)
        self.label_14.setObjectName(u"label_14")
        self.label_14.setFont(font1)

        self.horizontalLayout_9.addWidget(self.label_14)

        self.cbox_table_smp = QComboBox(self.tab_table)
        self.cbox_table_smp.setObjectName(u"cbox_table_smp")
        sizePolicy1.setHeightForWidth(self.cbox_table_smp.sizePolicy().hasHeightForWidth())
        self.cbox_table_smp.setSizePolicy(sizePolicy1)
        self.cbox_table_smp.setFont(font2)

        self.horizontalLayout_9.addWidget(self.cbox_table_smp)

        self.label_13 = QLabel(self.tab_table)
        self.label_13.setObjectName(u"label_13")
        self.label_13.setFont(font1)

        self.horizontalLayout_9.addWidget(self.label_13)

        self.cbox_table_chr = QComboBox(self.tab_table)
        self.cbox_table_chr.setObjectName(u"cbox_table_chr")
        sizePolicy1.setHeightForWidth(self.cbox_table_chr.sizePolicy().hasHeightForWidth())
        self.cbox_table_chr.setSizePolicy(sizePolicy1)
        self.cbox_table_chr.setFont(font2)

        self.horizontalLayout_9.addWidget(self.cbox_table_chr)

        self.horizontalSpacer_3 = QSpacerItem(118, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_9.addItem(self.horizontalSpacer_3)

        self.btn_preview_table = QPushButton(self.tab_table)
        self.btn_preview_table.setObjectName(u"btn_preview_table")
        self.btn_preview_table.setFont(font2)

        self.horizontalLayout_9.addWidget(self.btn_preview_table)

        self.btn_export_table = QPushButton(self.tab_table)
        self.btn_export_table.setObjectName(u"btn_export_table")
        self.btn_export_table.setFont(font2)

        self.horizontalLayout_9.addWidget(self.btn_export_table)

        self.verticalLayout_4.addLayout(self.horizontalLayout_9)

        self.gridLayout_4.addLayout(self.verticalLayout_4, 0, 0, 1, 1)

        self.tabOutput.addTab(self.tab_table, "")

        self.gridLayout_5.addWidget(self.tabOutput, 0, 0, 1, 1)

        self.gridLayout_6.addWidget(self.grpOutput, 4, 0, 1, 1)

        self.grpInput = QGroupBox(frmMain)
        self.grpInput.setObjectName(u"grpInput")
        self.grpInput.setFont(font)
        self.gridLayout = QGridLayout(self.grpInput)
        self.gridLayout.setObjectName(u"gridLayout")
        self.verticalLayout = QVBoxLayout()
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.horizontalLayout = QHBoxLayout()
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.label = QLabel(self.grpInput)
        self.label.setObjectName(u"label")
        sizePolicy2 = QSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Preferred)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy2)
        self.label.setFont(font1)

        self.horizontalLayout.addWidget(self.label)

        self.text_genome_file_path = DragLineEdit(self.grpInput)
        self.text_genome_file_path.setObjectName(u"text_genome_file_path")
        self.text_genome_file_path.setFont(font2)

        self.horizontalLayout.addWidget(self.text_genome_file_path)

        self.btn_load_genome = QPushButton(self.grpInput)
        self.btn_load_genome.setObjectName(u"btn_load_genome")
        sizePolicy.setHeightForWidth(self.btn_load_genome.sizePolicy().hasHeightForWidth())
        self.btn_load_genome.setSizePolicy(sizePolicy)
        self.btn_load_genome.setFont(font)

        self.horizontalLayout.addWidget(self.btn_load_genome)

        self.verticalLayout.addLayout(self.horizontalLayout)

        self.horizontalLayout_2 = QHBoxLayout()
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.label_3 = QLabel(self.grpInput)
        self.label_3.setObjectName(u"label_3")
        self.label_3.setFont(font1)

        self.horizontalLayout_2.addWidget(self.label_3)

        self.text_seq_depth_file_path = DragLineEdit(self.grpInput)
        self.text_seq_depth_file_path.setObjectName(u"text_seq_depth_file_path")
        self.text_seq_depth_file_path.setFont(font2)

        self.horizontalLayout_2.addWidget(self.text_seq_depth_file_path)

        self.btn_load_seqdepth = QPushButton(self.grpInput)
        self.btn_load_seqdepth.setObjectName(u"btn_load_seqdepth")
        self.btn_load_seqdepth.setFont(font)

        self.horizontalLayout_2.addWidget(self.btn_load_seqdepth)

        self.verticalLayout.addLayout(self.horizontalLayout_2)

        self.horizontalLayout_3 = QHBoxLayout()
        self.horizontalLayout_3.setObjectName(u"horizontalLayout_3")
        self.label_2 = QLabel(self.grpInput)
        self.label_2.setObjectName(u"label_2")
        self.label_2.setFont(font1)

        self.horizontalLayout_3.addWidget(self.label_2)

        self.text_mosdepth_folder_path = DragLineEdit(self.grpInput)
        self.text_mosdepth_folder_path.setObjectName(u"text_mosdepth_folder_path")
        self.text_mosdepth_folder_path.setFont(font2)

        self.horizontalLayout_3.addWidget(self.text_mosdepth_folder_path)

        self.btn_load_mosdepth = QPushButton(self.grpInput)
        self.btn_load_mosdepth.setObjectName(u"btn_load_mosdepth")
        self.btn_load_mosdepth.setFont(font)

        self.horizontalLayout_3.addWidget(self.btn_load_mosdepth)

        self.verticalLayout.addLayout(self.horizontalLayout_3)

        self.horizontalLayout_4 = QHBoxLayout()
        self.horizontalLayout_4.setObjectName(u"horizontalLayout_4")
        self.label_8 = QLabel(self.grpInput)
        self.label_8.setObjectName(u"label_8")
        self.label_8.setFont(font1)

        self.horizontalLayout_4.addWidget(self.label_8)

        self.text_smp_list_file_path = DragLineEdit(self.grpInput)
        self.text_smp_list_file_path.setObjectName(u"text_smp_list_file_path")
        self.text_smp_list_file_path.setFont(font2)

        self.horizontalLayout_4.addWidget(self.text_smp_list_file_path)

        self.btn_load_smp_list = QPushButton(self.grpInput)
        self.btn_load_smp_list.setObjectName(u"btn_load_smp_list")
        self.btn_load_smp_list.setFont(font)

        self.horizontalLayout_4.addWidget(self.btn_load_smp_list)

        self.verticalLayout.addLayout(self.horizontalLayout_4)

        self.horizontalLayout_5 = QHBoxLayout()
        self.horizontalLayout_5.setObjectName(u"horizontalLayout_5")
        self.label_4 = QLabel(self.grpInput)
        self.label_4.setObjectName(u"label_4")
        self.label_4.setFont(font1)

        self.horizontalLayout_5.addWidget(self.label_4)

        self.text_gene_list_file_path = DragLineEdit(self.grpInput)
        self.text_gene_list_file_path.setObjectName(u"text_gene_list_file_path")
        self.text_gene_list_file_path.setFont(font2)

        self.horizontalLayout_5.addWidget(self.text_gene_list_file_path)

        self.btn_load_gene_list = QPushButton(self.grpInput)
        self.btn_load_gene_list.setObjectName(u"btn_load_gene_list")
        self.btn_load_gene_list.setFont(font)

        self.horizontalLayout_5.addWidget(self.btn_load_gene_list)

        self.verticalLayout.addLayout(self.horizontalLayout_5)

        self.gridLayout.addLayout(self.verticalLayout, 0, 0, 1, 1)

        self.gridLayout_6.addWidget(self.grpInput, 0, 0, 1, 1)

        self.btn_run = QPushButton(frmMain)
        self.btn_run.setObjectName(u"btn_run")
        font3 = QFont()
        font3.setFamilies([u"JetBrains Mono"])
        font3.setPointSize(16)
        font3.setBold(True)
        self.btn_run.setFont(font3)

        self.gridLayout_6.addWidget(self.btn_run, 3, 0, 1, 1)

        self.btn_check = QPushButton(frmMain)
        self.btn_check.setObjectName(u"btn_check")
        self.btn_check.setFont(font3)

        self.gridLayout_6.addWidget(self.btn_check, 1, 0, 1, 1)

        self.gridLayout_7.addLayout(self.gridLayout_6, 0, 0, 1, 1)

        self.retranslateUi(frmMain)

        self.tabOutput.setCurrentIndex(0)

        QMetaObject.connectSlotsByName(frmMain)

    # setupUi

    def retranslateUi(self, frmMain):
        frmMain.setWindowTitle(QCoreApplication.translate("frmMain", u"popCNV", None))
        self.grpSettings.setTitle(QCoreApplication.translate("frmMain", u"Settings", None))
        self.label_5.setText(QCoreApplication.translate("frmMain", u"Wild group:", None))
        self.label_6.setText(QCoreApplication.translate("frmMain", u"Threads:", None))
        self.label_7.setText(QCoreApplication.translate("frmMain", u"Work path: ", None))
        self.text_wrkdir.setText("")
        self.btn_load_wrkdir.setText(QCoreApplication.translate("frmMain", u"...", None))
        self.grpOutput.setTitle(QCoreApplication.translate("frmMain", u"Output", None))
        self.label_9.setText(QCoreApplication.translate("frmMain", u"Type:", None))
        self.label_10.setText(QCoreApplication.translate("frmMain", u"Sample:", None))
        self.label_11.setText(QCoreApplication.translate("frmMain", u"Chr:", None))
        self.btn_draw_pic.setText(QCoreApplication.translate("frmMain", u"DrawPic", None))
        self.btn_export_pic.setText(QCoreApplication.translate("frmMain", u"Export", None))
        self.tabOutput.setTabText(self.tabOutput.indexOf(self.tab_pic),
                                  QCoreApplication.translate("frmMain", u"Pictures", None))
        self.label_12.setText(QCoreApplication.translate("frmMain", u"Type:", None))
        self.label_14.setText(QCoreApplication.translate("frmMain", u"Sample:", None))
        self.label_13.setText(QCoreApplication.translate("frmMain", u"Chr:", None))
        self.btn_preview_table.setText(QCoreApplication.translate("frmMain", u"Preview", None))
        self.btn_export_table.setText(QCoreApplication.translate("frmMain", u"Export", None))
        self.tabOutput.setTabText(self.tabOutput.indexOf(self.tab_table),
                                  QCoreApplication.translate("frmMain", u"Tables", None))
        self.grpInput.setTitle(QCoreApplication.translate("frmMain", u"Input files", None))
        self.label.setText(QCoreApplication.translate("frmMain", u"Genome:     ", None))
        self.text_genome_file_path.setText("")
        self.btn_load_genome.setText(QCoreApplication.translate("frmMain", u"...", None))
        self.label_3.setText(QCoreApplication.translate("frmMain", u"Seq depth:  ", None))
        self.text_seq_depth_file_path.setText("")
        self.btn_load_seqdepth.setText(QCoreApplication.translate("frmMain", u"...", None))
        self.label_2.setText(QCoreApplication.translate("frmMain", u"Mosdepth:   ", None))
        self.text_mosdepth_folder_path.setText("")
        self.btn_load_mosdepth.setText(QCoreApplication.translate("frmMain", u"...", None))
        self.label_8.setText(QCoreApplication.translate("frmMain", u"Sample list:", None))
        self.text_smp_list_file_path.setText("")
        self.btn_load_smp_list.setText(QCoreApplication.translate("frmMain", u"...", None))
        self.label_4.setText(QCoreApplication.translate("frmMain", u"Gene list:  ", None))
        self.text_gene_list_file_path.setText("")
        self.btn_load_gene_list.setText(QCoreApplication.translate("frmMain", u"...", None))
        self.btn_run.setText(QCoreApplication.translate("frmMain", u"Run", None))
        self.btn_check.setText(QCoreApplication.translate("frmMain", u"Check", None))
    # retranslateUi
