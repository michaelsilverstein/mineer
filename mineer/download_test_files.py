"""
Download test files
"""

import os, requests

def download(outdir=None):
    # URLs to test files
    root = 'https://raw.githubusercontent.com/michaelsilverstein/database/main/mineer/'
    files = ['SRR9660307_1.fastq', 'SRR9660307_2.fastq', 'SRR9660321_1.fastq', 'SRR9660321_2.fastq']
    urls = list(map(lambda f: os.path.join(root, f), files))

    if not outdir:
        outdir = 'sample_files'
    
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
        
    for url in urls:
        filename = os.path.basename(url)
        outpath = os.path.join(outdir, filename)
        print(f'Downloading {filename} to {outpath}')
        with requests.get(url) as resp:
            with open(outpath, 'wb') as fh:
                fh.write(resp.content)