from PySide6.QtCore import QAbstractTableModel, Qt


class TableModel(QAbstractTableModel):
    def __init__(self, data, header_data):
        super(TableModel, self).__init__()
        self.__data = data
        self.__header_data = header_data

    def data(self, index, role):
        if role == Qt.DisplayRole:
            return self.__data[index.row()][index.column()]

    def load_data(self, data, header_data):
        self.__data = data
        self.__header_data = header_data

    def rowCount(self, index):
        return len(self.__data)

    def columnCount(self, index):
        return len(self.__data[0])

    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.__header_data[section]
