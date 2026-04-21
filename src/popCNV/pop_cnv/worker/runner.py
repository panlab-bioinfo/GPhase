from os import popen
from pop_cnv.io.message import Message


class Runner:
    def __init__(self):
        self.__cmd = ""
        self.__res = ""

    def set_cmd(self, cmd):
        self.__cmd = cmd

    def run(self):
        try:
            with popen(self.__cmd, 'r') as fin:
                res = []
                for _ in fin:
                    res.append(_)
            self.__res = '\n'.join(res)
            return 1
        except Exception as e:
            Message.error(repr(e))

    def get_result(self):
        return self.__res
