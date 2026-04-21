#!/usr/bin/env python3
import sys
from os import path
from pathos import multiprocessing
from PySide6.QtWidgets import QApplication
from PySide6.QtGui import QIcon
from PySide6.QtCore import QCoreApplication, Qt
from pop_cnv.pipeline.gui_pipeline import GPopCNV

if __name__ == "__main__":
    multiprocessing.freeze_support()
    QCoreApplication.setAttribute(Qt.AA_ShareOpenGLContexts)
    app = QApplication([])
    icon_path = path.join(path.dirname(path.abspath(__file__)), "resources/icon.png")
    app.setWindowIcon(QIcon(icon_path))
    main_window = GPopCNV()
    main_window.show()
    sys.exit(app.exec())
