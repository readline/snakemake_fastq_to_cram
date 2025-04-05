from scripts.utils import allocated,ignore

rule cramqc__flagstat:
    input:
        cram= join(config['workdir'], "02.Alignment", "Level3", "{sample}", "{sample}.BQSR.cram"),
    output:
        flagstat = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Flagstat", "Flagstat.json"),
    log: 
        out = join(config['pipelinedir'], "logs", "cramqc__flagstat", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cramqc__flagstat", "{sample}.e"),
    params:
        folder = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Flagstat"),
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
        cd {params.folder}
        samtools flagstat -@ {resources.cpus_per_task} -O json {input.cram} > {output.flagstat} 2>{log.err}
        """

rule cramqc__collect_quality_yield_metrics:
    input:
        cram = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "{sample}.BQSR.cram"),
    output:
        metrix = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Metrics", "QualityYield.metrics"),
    log:
        out = join(config['pipelinedir'], "logs", "cramqc__collect_quality_yield_metrics", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cramqc__collect_quality_yield_metrics", "{sample}.e"),
    resources:
        cpus_per_task = 4,
        mem = '4G',
        runtime = '7d',
        partition = 'defq',
    conda:
        config['conda']['align']
    # container:
        # config['container']['gatk']
    shell:
        "gatk --java-options \"-Xms2000m -Xmx3000m\" "
        "  CollectQualityYieldMetrics"
        "  -I {input.cram}"
        "  -O {output.metrix}"
        "  --CREATE_MD5_FILE true"
        "  -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  > {log.out} 2> {log.err}\n"

rule cramqc__collect_wgs_metrics:
    input:
        cram = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "{sample}.BQSR.cram"),
    output:
        metrix   = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Metrics", "WGS.metrics"),
    log:
        out = join(config['pipelinedir'], "logs", "cramqc__collect_wgs_metrics", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cramqc__collect_wgs_metrics", "{sample}.e"),
    resources:
        cpus_per_task = 4,
        mem = '4G',
        runtime = '7d',
        partition = 'defq',
    conda:
        config['conda']['align']
    # container:
        # config['container']['gatk']
    shell:
        "gatk --java-options \"-Xms2000m -Xmx3000m\" "
        " CollectWgsMetrics"
        "  -I {input.cram}"
        "  -O {output.metrix}"
        "  -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  --INCLUDE_BQ_HISTOGRAM true"
        "  --INTERVALS {config[references][gatkbundle]}/wgs_coverage_regions.hg38.interval_list"
        "  --VALIDATION_STRINGENCY SILENT"
        "  --USE_FAST_ALGORITHM true"
        "  --CREATE_MD5_FILE true"
        "  --READ_LENGTH 151"
        "  > {log.out} 2> {log.err}\n"

rule cramqc__collect_all_reads_multiple_metrics:
    input:
        cram = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "{sample}.BQSR.cram"),
    output:
        status = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Metrics", "AllReadsMultiple.ok"),
    params:
        prefix = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Metrics", "AllReadsMultiple.metrics"),
    log:
        out = join(config['pipelinedir'], "logs", "cramqc__collect_all_reads_multiple_metrics", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cramqc__collect_all_reads_multiple_metrics", "{sample}.e"),
    resources:
        cpus_per_task = 4,
        mem = '8G',
        runtime = '7d',
        partition = 'defq',
    conda:
        config['conda']['align']
    # container:
        # config['container']['gatk']
    shell:
        "gatk --java-options \"-Xms5000m -Xmx6500m\" "
        "  CollectMultipleMetrics"
        "  -I {input.cram}"
        "  -O {params.prefix}"
        "  --ASSUME_SORTED true"
        "  --CREATE_MD5_FILE true"
        "  -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  --VALIDATION_STRINGENCY LENIENT"
        "  --PROGRAM null "
        "  --PROGRAM CollectBaseDistributionByCycle "
        "  --PROGRAM CollectInsertSizeMetrics "
        "  --PROGRAM MeanQualityByCycle "
        "  --PROGRAM QualityScoreDistribution "
        "  --METRIC_ACCUMULATION_LEVEL null "
        "  --METRIC_ACCUMULATION_LEVEL ALL_READS "
        " > {log.out} 2> {log.err}\n"
        "touch {output.status}"
        " >> {log.out} 2>> {log.err}\n"

rule cramqc__collect_read_groups_multiple_metrics:
    input:
        cram = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "{sample}.BQSR.cram"),
    output:
        status = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Metrics", "ReadGroupsMultiple.ok"),
    params:
        prefix = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Metrics", "ReadGroupsMultiple.metrics"),
    log:
        out = join(config['pipelinedir'], "logs", "cramqc__collect_read_groups_multiple_metrics", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cramqc__collect_read_groups_multiple_metrics", "{sample}.e"),
    resources:
        cpus_per_task = 4,
        mem = '8G',
        runtime = '7d',
        partition = 'defq',
    conda:
        config['conda']['align']
    # container:
        # config['container']['gatk']
    shell:
        "gatk --java-options \"-Xms5000m -Xmx6500m\" "
        "  CollectMultipleMetrics"
        "  -I {input.cram}"
        "  -O {params.prefix}"
        "  --ASSUME_SORTED true"
        "  --CREATE_MD5_FILE true"
        "  -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  --VALIDATION_STRINGENCY LENIENT"
        "  --PROGRAM null "
        "  --PROGRAM CollectBaseDistributionByCycle "
        "  --PROGRAM CollectInsertSizeMetrics "
        "  --PROGRAM CollectAlignmentSummaryMetrics "
        "  --PROGRAM QualityScoreDistribution "
        "  --METRIC_ACCUMULATION_LEVEL ALL_READS "
        "  --METRIC_ACCUMULATION_LEVEL READ_GROUP "
        " > {log.out} 2> {log.err}\n"
        "touch {output.status}"
        "  >> {log.out} 2>> {log.err}\n"

rule cramqc__collect_aggregation_metrics:
    input:
        cram = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "{sample}.BQSR.cram"),
    output:
        status   = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Metrics", "SM_LB_Aggregation.ok"),
    params:
        prefix = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Metrics", "SM_LB_Aggregation.metrics"),
    log:
        out = join(config['pipelinedir'], "logs", "cramqc__collect_aggregation_metrics", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cramqc__collect_aggregation_metrics", "{sample}.e"),
    resources:
        cpus_per_task = 4,
        mem = '8G',
        runtime = '7d',
        partition = 'defq',
    conda:
        config['conda']['align']
    # container:
        # config['container']['gatk']
    shell:
        "gatk --java-options \"-Xms5000m -Xmx6500m\" "
        "  CollectMultipleMetrics"
        "  -I {input.cram}"
        "  -O {params.prefix}"
        "  --ASSUME_SORTED true"
        "  --CREATE_MD5_FILE true"
        "  -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  --PROGRAM null "
        "  --PROGRAM CollectAlignmentSummaryMetrics "
        "  --PROGRAM CollectInsertSizeMetrics "
        "  --PROGRAM CollectSequencingArtifactMetrics "
        "  --PROGRAM QualityScoreDistribution "
        "  --METRIC_ACCUMULATION_LEVEL null "
        "  --METRIC_ACCUMULATION_LEVEL SAMPLE "
        "  --METRIC_ACCUMULATION_LEVEL LIBRARY "
        " >> {log.out} 2>> {log.err}\n"
        "touch {output.status}"
        " >> {log.out} 2>> {log.err}\n"

rule cramqc__mosdepth:
    input:
        cram = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "{sample}.BQSR.cram"),
    output:
        status = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Mosdepth", "{sample}.mosdepth.summary.txt"),
    params:
        prefix = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Mosdepth", "{sample}"),
    log:
        out = join(config['pipelinedir'], "logs", "cramqc__mosdepth", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cramqc__mosdepth", "{sample}.e"),
    resources:
        cpus_per_task = 16,
        mem = '32G',
        runtime = '7d',
        partition = 'defq',
    conda:
        config['conda']['align']
    # container:
        # config['container']['mosdepth']
    shell:
        "export MOSDEPTH_Q0=NO_COVERAGE \n"
        "export MOSDEPTH_Q1=LOW_COVERAGE \n"
        "export MOSDEPTH_Q2=CALLABLE \n"
        "export MOSDEPTH_Q3=HIGH_COVERAGE \n"
        "export MOSDEPTH_Q4=HIGH_COVERAGE_300 \n"
        "export MOSDEPTH_Q5=HIGH_COVERAGE_1000 \n"
        "mosdepth"
        "    --threads {resources.cpus_per_task}"
        "    --fasta {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "    --quantize 0:1:5:150:300:1000:"
        "    {params.prefix}"
        "    {input.cram}"
        " >> {log.out} 2>> {log.err}\n"

rule cramqc__fingerprint:
    input:
        cram = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "{sample}.BQSR.cram"),
    output:
        fingerprint = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Fingerprint", "fingerprint.vcf"),
        metrics     = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Fingerprint", "fingerprint.rg.metrics"),
    log:
        out = join(config['pipelinedir'], "logs", "cramqc__fingerprint", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cramqc__fingerprint", "{sample}.e"),
    resources:
        cpus_per_task = 4,
        mem = '8G',
        runtime = '7d',
        partition = 'defq',
    conda:
        config['conda']['align']
    # container:
        # config['container']['mosdepth']
    shell:
        "gatk --java-options \"-Xms5000m -Xmx6500m\" "
        "    ExtractFingerprint"
        "    -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "    -H {config[references][gatkbundle]}/Homo_sapiens_assembly38.haplotype_database.txt"
        "    -I {input.cram}"
        "    -O {output.fingerprint}"
        " > {log.out} 2> {log.err}\n"
        "gatk --java-options \"-Xms5000m -Xmx6500m\" "
        "    CrosscheckFingerprints"
        "    -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "    -H {config[references][gatkbundle]}/Homo_sapiens_assembly38.haplotype_database.txt"
        "    -I {input.cram}"
        "    -O {output.metrics}"
        " >> {log.out} 2>> {log.err}\n"

rule cramqc__summary_metrics:
    input:
        join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Flagstat", "Flagstat.json"),
        join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Metrics", "QualityYield.metrics"),
        join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Metrics", "WGS.metrics"),
        join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Metrics", "AllReadsMultiple.ok"),
        join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Metrics", "ReadGroupsMultiple.ok"),
        join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Metrics", "SM_LB_Aggregation.ok"),
        join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Mosdepth", "{sample}.mosdepth.summary.txt"),
        join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Fingerprint", "fingerprint.rg.metrics"),
    output:
        status   = join(config['workdir'], "02.Alignment", "Level3", "{sample}", "QC", "Summary.ok"),
    log:
        out = join(config['pipelinedir'], "logs", "cramqc__summary_metrics", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cramqc__summary_metrics", "{sample}.e"),
    resources:
        cpus_per_task = 1,
        mem = '1G',
        runtime = '2d',
        partition = 'defq',
    shell:
        "touch {output.status}"
