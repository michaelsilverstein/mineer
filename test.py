import multiprocessing as mp
from mineer.utils import File
from mineer.mineer import minEER
import random
from time import time

def mpRunMineer(read):
    trimpos_mineer = minEER(read.untrimmed.ee, read.mal, read.mae)
    return trimpos_mineer


if __name__ == '__main__':
    file = File('test_files/SRR9660307_1.fastq', 'f', 100, 1e-2)
    reads = file.reads
    for r in reads:
        r.untrimmed

    start = time()
    for read in reads:
        read.runMineer()
    print('Single thread')
    print(f'{time() - start:.3e} seconds')
    print([r.trimpos_mineer for r in reads[:10]])

    start = time()
    pool = mp.Pool()
    print(pool)
    jobs = pool.map_async(mpRunMineer, reads, 1000).get()
    print('Multiprocessed')
    print(f'{time() - start:.3e} seconds')
    print(jobs[:10])

    # nreads = 5000
    # random.shuffle(reads)

    # read_count = 0
    # passing_reads =[]
    
    # i = 0
    # pool = mp.Pool()
    
    # chunk = reads[(i * 1000) : (i + 1) * 1000]
    # pool.map_async(mpRunMineer, chunk, 50)
    # for read in chunk:
    #     if read.pass_qc_mineer:
    #         passing_reads.append(read)
    # print(len(passing_reads))