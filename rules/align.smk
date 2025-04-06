from scripts.utils import allocated,ignore

itv4 = ['%.4d'%(itv) for itv in range(1, int(config['intervals'])+1)]

rule Bwa_mem:
    """
    Align reads to reference genome using BWA-MEM.
    """
    input:
        read1 = join(config['workdir'], "01.clean_data", "{run}", "{run}.R1.fastq.gz"),
        read2 = join(config['workdir'], "01.clean_data", "{run}", "{run}.R2.fastq.gz"),
    output:
        bam = temp(join(config['workdir'], "02.Alignment", "Level1", "{run}", "{run}.sort.bam")),
        bai = temp(join(config['workdir'], "02.Alignment", "Level1", "{run}", "{run}.sort.bai")),
    log:
        out = join(config['pipelinedir'], "logs", "Bwa_mem", "{run}.o"),
        err = join(config['pipelinedir'], "logs", "Bwa_mem", "{run}.e"),
    params:
        bwa=lambda wildcards:' -K 10000000 -R "@RG\\tID:{}\\tLB:{}\\tPL:illumina\\tPU:{}\\tSM:{}" '.format(
                                wildcards.run, rundic[wildcards.run]['Lib'],  wildcards.run,  rundic[wildcards.run]['Sample']),
        tmpdir=lambda wildcards: join(config['workdir'], "02.Alignment", "Level1", wildcards.run, wildcards.run + ".tmp"),
    resources:
        cpus_per_task = 24,
        mem = '60G',
        runtime = '2d',
        partition = 'defq',
    conda:
        config['conda']['align']
    # container:
        # config['container']['gatk']
    shell:
        "mkdir -p {params.tmpdir} \n"
        "bwa mem"
        "  -t {resources.cpus_per_task}"
        "  {params.bwa}"
        "  {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  {input.read1}"
        "  {input.read2}"
        "  2> {log.err} |"
        "gatk SortSam --java-options \"-Xms30G -Xmx30G -XX:ParallelGCThreads=2 -Djava.io.tmpdir={params.tmpdir}\""
        "  --MAX_RECORDS_IN_RAM 5000000"
        "  -I /dev/stdin"
        "  --SORT_ORDER coordinate"
        "  -O /dev/stdout"
        "  --TMP_DIR {params.tmpdir}"
        "  2>> {log.err} |"
        "gatk SetNmMdAndUqTags --java-options \" -Xms30G -Xmx30G -XX:ParallelGCThreads=2 -Djava.io.tmpdir={params.tmpdir}\" "
        "  --INPUT /dev/stdin "
        "  --OUTPUT {output.bam} "
        "  --REFERENCE_SEQUENCE {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  --CREATE_INDEX true"
        "  --TMP_DIR {params.tmpdir}"
        "  >> {log.out} 2>> {log.err} \n"
        "rm -rf {params.tmpdir}"

rule Markdup:
    input:
        bam = lambda wildcards: expand(join(config['workdir'], "02.Alignment", "Level1", "{run}", "{run}.sort.bam"), run=libdic[wildcards.lib]),
        bai = lambda wildcards: expand(join(config['workdir'], "02.Alignment", "Level1", "{run}", "{run}.sort.bai"), run=libdic[wildcards.lib]),
    output:
        bam    =temp(join(config['workdir'], "02.Alignment", "Level2", "{lib}", "{lib}.sort.md.bam")),
        bai    =temp(join(config['workdir'], "02.Alignment", "Level2", "{lib}", "{lib}.sort.md.bam.bai")),
        metrics=join(config['workdir'], "02.Alignment", "Level2", "{lib}", "{lib}.sort.md.metrics"),
    log:
        out = join(config['pipelinedir'], "logs", "Markdup", "{lib}.o"),
        err = join(config['pipelinedir'], "logs", "Markdup", "{lib}.e"),
    params:
        tmpdir=lambda wildcards: join(config['workdir'], "02.Alignment", "Level2", wildcards.lib, wildcards.lib + ".tmp"),
    resources:
        cpus_per_task = 8,
        mem = '60G',
        runtime = '2d',
        partition = 'defq',
    conda:
        config['conda']['align']
    # container:
        # config['container']['gatk']
    shell:
        "mkdir -p {params.tmpdir} \n"
        "gatk --java-options \"-Xms30G -Xmx30G -XX:ParallelGCThreads=2 -Djava.io.tmpdir={params.tmpdir}\""
        "  MarkDuplicates"
        "  $(for f in {input.bam}; do echo -n \" -I $f\"; done)"
        "  -O {output.bam}"
        "  -M {output.metrics}"
        "  --MAX_RECORDS_IN_RAM 2000000"
        "  --VALIDATION_STRINGENCY SILENT"
        "  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500"
        "  --TMP_DIR {params.tmpdir}"
        "  > {log.out} 2> {log.err} \n"
        "sambamba index -t {resources.cpus_per_task} {output.bam} >> {log.out} 2>> {log.err} \n"
        "rm -rf {params.tmpdir}"
        
rule Merge_level3:
    input:
        bam = lambda wildcards: expand(join(config['workdir'], "02.Alignment", "Level2", "{lib}", "{lib}.sort.md.bam"), lib=sampledic[wildcards.sample]),
        bai = lambda wildcards: expand(join(config['workdir'], "02.Alignment", "Level2", "{lib}", "{lib}.sort.md.bam.bai"), lib=sampledic[wildcards.sample]),
    output:
        bam=temp(join(config['workdir'], "02.Alignment", "Level3", "{sample}", "{sample}.sort.md.bam")),
        bai=temp(join(config['workdir'], "02.Alignment", "Level3", "{sample}", "{sample}.sort.md.bam.bai")),
    log:
        out = join(config['pipelinedir'], "logs", "Merge_level3", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "Merge_level3", "{sample}.e"),
    resources:
        cpus_per_task = 8,
        mem = '16G',
        runtime = '2d',
        partition = 'defq',
    conda:
        config['conda']['align']
    # container:
        # config['container']['gatk']
    shell:
        """
        if [ $(echo {input.bam} | wc -w) -gt 1 ]; then
            samtools merge \
                -@ {resources.cpus_per_task} \
                {output.bam} \
                {input.bam} > {log.out} 2> {log.err}
            samtools index \
                -@ {resources.cpus_per_task} \
                {output.bam} >> {log.out} 2>> {log.err}
        else
            mv {input.bam} {output.bam} > {log.out} 2> {log.err}
            samtools index \
                -@ {resources.cpus_per_task} \
                {output.bam} >> {log.out} 2>> {log.err}
        fi
        """

rule BQSR_ht:
    input:
        bam = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "{sample}.sort.md.bam"),
        bai = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "{sample}.sort.md.bam.bai"),
    output:
        metrics=temp(join(config['workdir'], "02.Alignment", "Level3", "{sample}", "{sample}.BQSR.metrics")),
        bam    =temp(join(config['workdir'], "02.Alignment", "Level3", "{sample}", "{sample}.BQSR.bam")),
    params:
        itvs = " ".join(itv4),
        tmpdir=lambda wildcards: join(config['workdir'], "02.Alignment", "Level3", wildcards.sample, wildcards.sample + ".tmp"),
    log:
        out = join(config['pipelinedir'], "logs", "BQSR", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "BQSR", "{sample}.e"),
    resources:
        cpus_per_task = 32,
        mem = '128G',
        runtime = '2d',
        partition = 'defq',
    conda:
        config['conda']['align']
    # container:
        # config['container']['gatk']
    shell:
        """mkdir -p {params.tmpdir} \
            > {log.out} 2> {log.err}
        echo "{params.itvs}" | tr ' ' '\\n' | parallel -j 16 '
            chunk={{}}
            mkdir -p {params.tmpdir}/$chunk \
                >> {log.out} 2>> {log.err} 
            gatk --java-options \"-Xms4G -Xmx4G -XX:ParallelGCThreads=2 -Djava.io.tmpdir={params.tmpdir}/$chunk\" \
                BaseRecalibrator \
                -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta \
                -I {input.bam} \
                -O {params.tmpdir}/$chunk/$chunk.BQSR.metrics \
                --use-original-qualities \
                --known-sites {config[references][gatkbundle]}/Homo_sapiens_assembly38.dbsnp138.vcf \
                --known-sites {config[references][gatkbundle]}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
                --known-sites {config[references][gatkbundle]}/Homo_sapiens_assembly38.known_indels.vcf.gz \
                --intervals {config[references][gatkbundle]}/scattered_calling_intervals/temp_$chunk\_of_50.bed \
                >> {log.out} 2>> {log.err}
            gatk --java-options \"-Xms4G -Xmx4G -XX:ParallelGCThreads=2 -Djava.io.tmpdir={params.tmpdir}/$chunk\" \
                ApplyBQSR \
                --add-output-sam-program-record \
                -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta \
                -I {input.bam} \
                -O {params.tmpdir}/$chunk/$chunk.BQSR.bam \
                --use-original-qualities \
                -bqsr {params.tmpdir}/$chunk/$chunk.BQSR.metrics \
                --static-quantized-quals 10 \
                --static-quantized-quals 20 \
                --static-quantized-quals 30 \
                -L {config[references][gatkbundle]}/scattered_calling_intervals/temp_$chunk\_of_50.bed \
                >> {log.out} 2>> {log.err}'
        inputs=$( ls {params.tmpdir}/*/*.BQSR.metrics | sed 's/^/-I\ /g' )
        gatk --java-options \"-Xms3000m\" \
            GatherBQSRReports \
            $inputs \
            -O {output.metrics} \
            >> {log.out} 2>> {log.err}
        inputs=$( ls {params.tmpdir}/*/*.BQSR.bam | sed 's/^/-I\ /g' )
        gatk MergeSamFiles \
            $inputs \
            -O {output.bam} \
            --USE_THREADING true \
            >> {log.out} 2>> {log.err}
        """

rule Bam_to_cram:
    input:
        bam = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "{sample}.BQSR.bam"),
    output:
        cram= join(config['workdir'], "02.Alignment", "Level3", "{sample}", "{sample}.BQSR.cram"),
    log:
        out = join(config['pipelinedir'], "logs", "Bam_to_cram", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "Bam_to_cram", "{sample}.e"),
    resources:
        cpus_per_task = 8,
        mem = '32G',
        runtime = '2d',
        partition = 'normal_q',
    conda:
        config['conda']['align']
    # container:
        # config['container']['gatk']
    shell:
        "samtools view"
        "  -@ {resources.cpus_per_task}"
        "  --reference {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  --output-fmt CRAM"
        "  -o {output.cram}"
        "  {input.bam}"
        "  > {log.out} 2> {log.err}\n"
        "samtools index"
        "  -@ {resources.cpus_per_task}"
        "  {output.cram}"
        "  >> {log.out} 2>> {log.err}"
