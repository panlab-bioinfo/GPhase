from os import popen


class DepCheck:
    def __init__(self):
        pass

    @staticmethod
    def check():
        res = []
        with popen("mosdepth -h", 'r') as fin:
            for line in fin:
                res.append(line.strip())

        if res and res[0].startswith("mosdepth"):
            return True
        return False
