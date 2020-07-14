from builtins import object
import multiprocessing
import os


class Pool(object):
    def __init__(self, processes):
        self.size = processes
        self.rank = 0
        self.master_pid = os.getpid()

    def is_master(self):
        return self.master_pid == os.getpid()

    def map(self, function, args):
        pool = multiprocessing.Pool(self.size)
        results = pool.map(function, args)
        pool.close()
        return results

    def close(self):
        pass

    def bcast(self, data):
        return self.data

    def gather(self, data):
        return self.data


    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass        