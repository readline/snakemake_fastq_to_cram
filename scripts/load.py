#!/usr/bin/env python
# author: kaiyu062@gmail.com
# date: 2025-01-06
# Function to load samplesheet and return sample, library, run dictionaries.
# Can be tested with:
# python load.py samplesheet.tsv
import sys
import pandas as pd

def samplesheet(sspath):
    df = pd.read_csv(sspath, sep='\t', usecols=range(8))
    idx = []
    for i in df.index:
        idx.append('%s.%s.%s.%s'%(df.loc[i,'Sample'], df.loc[i,'Lib'], df.loc[i,'Flowcell'], df.loc[i,'Lane']))
    df.index=idx
    sample = {}
    lib = {}
    run = {}
    for i in df.index:
        ss = str(df.loc[i,'Sample'])
        li = str(df.loc[i,'Lib'])
        if ss not in sample:
            sample[ss] = []
        if li not in lib:
            lib[li] = []
        if li not in sample[ss]:
            sample[ss].append(li)
        if i not in lib[li]:
            lib[li].append(i)
        run[i] = {n:df.loc[i,n] for n in df.columns}
    return sample,lib,run

if __name__ == '__main__':
    print(samplesheet(sys.argv[1]))