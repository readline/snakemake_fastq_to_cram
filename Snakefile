import os, yaml
import pandas as pd
from os.path import join
from scripts.load import samplesheet

snakedir = os.path.dirname(workflow.snakefile)
print('Snakefile dir: ', snakedir)
configfile: join(snakedir, 'config.yaml')
sampledic, libdic, rundic = samplesheet(config['samplesheet'])
workdir: config['workdir']

print(sampledic, libdic, rundic)

# Load cluster config
include: join(config['pipelinedir'], "rules", "prepare.smk")
include: join(config['pipelinedir'], "rules", "align.smk")
include: join(config['pipelinedir'], "rules", "stats.smk")

rule all:
    input:
        expand(join(config['workdir'], "02.Alignment", "Level3", "{sample}", "{sample}.BQSR.cram"),
            sample=sampledic),
        expand(join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Summary.ok"),
            sample=sampledic),
