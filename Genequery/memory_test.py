import myvariant
import time
import os
import psutil
from pyliftover import LiftOver
import math
import re
import csv
import datetime
import pandas as pd


def read_in_chunks(file_object, chunk_size=1024):
    """Lazy function (generator) to read a file piece by piece.
    Default chunk size: 1k."""
    while True:
        data = file_object.read(chunk_size)
        if not data:
            break
        yield data

def get_process_memory():
    pid = os.getpid()
    print(pid)
    ps = psutil.Process(pid)
    memoryUse = ps.memory_info()
    print('memoryUse.rss : ',memoryUse.rss)
    return memoryUse.rss

#for piece in read_in_chunks(f):
#    print(len(piece))

num = 0

mem_before = get_process_memory()
start = time.time()
f = open('23andme_large.txt')
elapsed_time = time.time() - start
mem_after = get_process_memory()
print(num)
print("{}: memory before: {:,}, after: {:,}, consumed: {:,}; exec time: {}".format(
            'test_func',
            mem_before, mem_after, mem_after - mem_before,
            elapsed_time))
