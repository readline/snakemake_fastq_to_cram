from scripts.utils import allocated,ignore

rule Fastp:
    """
    Quality-control step to assess sequencing quality and remove adapter sequences.
    FastP performs quality control, adapter trimming, and filtering of raw sequencing data.
    @Input:
        Raw FastQ files (paired-end)
    @Output:
        Trimmed FastQ files and QC reports
    """
    input:
        r1 = lambda wildcards: rundic[wildcards.run]['Read1'],
        r2 = lambda wildcards: rundic[wildcards.run]['Read2'],
    output:
        r1 = temp(join(config['workdir'], "01.clean_data", "{run}", "{run}.R1.fastq.gz")),
        r2 = temp(join(config['workdir'], "01.clean_data", "{run}", "{run}.R2.fastq.gz")),
        html = join(config['workdir'], "01.clean_data", "{run}", "{run}.fastp.html"),
        json = join(config['workdir'], "01.clean_data", "{run}", "{run}.fastp.json"),
    log:
        out = join(config['pipelinedir'], "logs", "Fastp", "{run}.o"),
        err = join(config['pipelinedir'], "logs", "Fastp", "{run}.e"),
    params:
        folder = directory(join(config['workdir'],"01.clean_data","{run}")),
    resources:
        cpus_per_task = 12,
        mem = '24G',
        slurm_partition = 'defq',
        runtime = '2d', 
    conda:
        config['conda']['align']
    # container:
        # config['container']['fastp']
    shell:
        "cd {params.folder} \n"
        "fastp"
        " -i {input.r1}" 
        " -I {input.r2}" 
        " -o {output.r1}" 
        " -O {output.r2}" 
        " -t 1"
        " -T 1"
        " --dont_eval_duplication"
        " -h {output.html}" 
        " -j {output.json}" 
        " -w {resources.cpus_per_task}"
        "  > {log.out} 2> {log.err}"
        
