import matplotlib

matplotlib.use("qt5agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
from numpy import log2, errstate


class PlotCanvas(FigureCanvas):

    def __init__(self, width=20, height=8, dpi=300):
        fig = plt.figure(figsize=(width, height), dpi=dpi, tight_layout=True)
        self.plt = plt
        super(PlotCanvas, self).__init__(fig)


class PlotContent:

    def __init__(self):
        self.figure_content = PlotCanvas()

    def gene_line_graph(self, data_db):
        self.figure_content.plt.clf()
        self.figure_content.plt.rcParams["figure.figsize"] = (20, 8)
        self.figure_content.plt.rcParams["figure.dpi"] = 300
        try:
            xticks = []
            xlabels = []
            is_first = True
            for smp in data_db:
                off_set = 0
                X = []
                Y = []
                data = []
                for chrn in sorted(data_db[smp]):
                    for idx in range(len(data_db[smp][chrn]["X"])):
                        x = data_db[smp][chrn]["X"][idx]
                        y = data_db[smp][chrn]["Y"][idx]
                        data.append([x + off_set, y])
                    if is_first:
                        xticks.append(off_set + max(data_db[smp][chrn]["X"]) / 2.)
                        xlabels.append(chrn)
                    off_set += max(data_db[smp][chrn]["X"])
                is_first = False
                X = [v[0] for v in sorted(data)]
                Y = [v[1] for v in sorted(data)]
                self.figure_content.plt.plot(X, Y, label=smp, lw=1)
            self.figure_content.plt.legend()
            self.figure_content.plt.xticks(xticks, xlabels, rotation=-45)
        except Exception as e:
            return str(e)
        return "Success"

    def gene_manhattan_graph(self, data_db):
        self.figure_content.plt.clf()
        try:
            colors = ['steelblue', 'darkorange']
            off_set = 0
            max_x = 0
            cidx = 0
            xticks = []
            xlabels = []
            for chrn in sorted(data_db):
                data = []
                for idx in range(len(data_db[chrn]["X"])):
                    x = data_db[chrn]["X"][idx]
                    y = data_db[chrn]["Y"][idx]
                    data.append([x + off_set, y])
                X = [v[0] for v in sorted(data)]
                Y = [v[1] for v in sorted(data)]
                with errstate(divide='ignore'):
                    self.figure_content.plt.scatter(X, -log2(Y), s=3, color=colors[cidx % 2])
                max_x = max(max_x, max(X))
                cidx += 1
                xticks.append(off_set + max(data_db[chrn]["X"]) / 2.)
                xlabels.append(chrn)
                off_set += max(data_db[chrn]["X"])
            self.figure_content.plt.plot([0, max_x], [-log2(.05), -log2(.05)],
                                         color='blue', linestyle=':', lw=1)
            self.figure_content.plt.plot([0, max_x], [-log2(.01), -log2(.01)],
                                         color='red', linestyle=':', lw=1)
            self.figure_content.plt.xticks(xticks, xlabels, rotation=-45)
            self.figure_content.plt.xlim(0, max_x)
        except Exception as e:
            return str(e)
        return "Success"
