import os
import sys
import functools


# GENERIC DATA

ORGANISMS = {
    'rn4': 'rat',
    'rn5': 'rat',
    'mm9': 'mouse',
    'mm10': 'mouse',
    'hg19': 'human',
    'hg38': 'human'
}

INCLUDE_CHRS = {
    'hg19': ['chr{}'.format(i) for i in range(1, 23)] + ["chrX", "chrY"],
    'hg38': ['chr{}'.format(i) for i in range(1, 23)] + ["chrX", "chrY"],
    'mm9': ['chr{}'.format(i) for i in range(1, 20)] + ["chrX", "chrY"],
    'mm10': ['chr{}'.format(i) for i in range(1, 20)] + ["chrX", "chrY"],
    'rn4': ['chr{}'.format(i) for i in range(1, 21)] + ["chrX", "chrY"],
    'rn5': ['chr{}'.format(i) for i in range(1, 21)] + ["chrX", "chrY"]
}

MACS2_GENOME_SIZE = {
    'rn4': 'mm',
    'rn5': 'mm',
    'mm9': 'mm',
    'mm10': 'mm',
    'hg19': 'hs',
    'hg38': 'hs'
}

# Helper functions

def get_genome(sample):
    return(config['sample_genome'][sample])


# RESULT PATHS

prefix_results = functools.partial(os.path.join, config['results_dir'])
CONCAT_READS_DIR = prefix_results('concat_reads')
ALIGN_DIR = prefix_results('aligned')
PRUNE_DIR = prefix_results('pruned')
BIGWIG_DIR = prefix_results('bigwig') #formerly DISP_DIR
BEDGRAPH_DIR = prefix_results('bedgraph')
MACS2_NARROW_DIR = prefix_results('macs2_narrowpeak')
MACS2_BROAD_DIR = prefix_results('macs2_broadpeak')
SEACR_DIR = prefix_results('seacr_peak')
HOMERTAG_DIR = prefix_results('homer_tag')
HOMERPEAK_DIR = prefix_results('homer_peak')
HOMERMOTIF_DIR = prefix_results('homer_motif')
MACS2_NARROW_HOMERMOTIF_DIR = prefix_results('macs2_narrowpeak_homer_motif')
MACS2_BROAD_HOMERMOTIF_DIR = prefix_results('macs2_broadpeak_homer_motif')
SEACR_STRINGENT_HOMERMOTIF_DIR = prefix_results('seacr_stringentpeak_homer_motif')
SEACR_INTERMEDIATE_HOMERMOTIF_DIR = prefix_results('seacr_intermediatepeak_homer_motif')
SEACR_RELAXED_HOMERMOTIF_DIR = prefix_results('seacr_relaxedpeak_homer_motif')
SCRIPTS_DIR = os.path.join(workflow.basedir, 'scripts')
CONDA_ENVS_DIR = prefix_results('conda_envs')

# Set workdir - Snakemake will be run from this location.
workdir:
    config['results_dir']

###
# Conditional inputs to all rule
if config.get('sample_input'):
    if config.get('sample_homer_fmg_genome'):
        all_input = [
            expand(os.path.join(BIGWIG_DIR, "{sample}.1m.bw"), sample=config['sample_paths'].keys()), #Create bigwigs for all samples
            expand(os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.BLfiltered.hpeaks"), paramset=config['homer_findPeaks_params'].keys(), sample=config['sample_input'].keys()), #Call peaks for all samples with matched inputs
            expand(os.path.join(HOMERMOTIF_DIR, "{paramset}", "{sample}"), paramset=config['homer_findPeaks_params'].keys(), sample=config['sample_homer_fmg_genome'].keys()), #Homermotifs for all samples with a specified genome for homer findMotifsGenome
            expand(os.path.join(MACS2_NARROW_HOMERMOTIF_DIR, "{paramset}", "{sample}"), paramset=config['homer_findPeaks_params'].keys(), sample=config['sample_homer_fmg_genome'].keys()), #Homermotifs for all samples with a specified genome for homer findMotifsGenome
            expand(os.path.join(MACS2_BROAD_HOMERMOTIF_DIR, "{paramset}", "{sample}"), paramset=config['homer_findPeaks_params'].keys(), sample=config['sample_homer_fmg_genome'].keys()), #Homermotifs for all samples with a specified genome for homer findMotifsGenome
            expand(os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.BLfiltered.narrowPeak"), paramset=config['macs2_broad_params'].keys(), sample=config['sample_paths'].keys()),
            expand(os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.summits.BLfiltered.bed"), paramset=config['macs2_broad_params'].keys(), sample=config['sample_paths'].keys()),
            expand(os.path.join(MACS2_NARROW_DIR, "{paramset}", "{sample}.BLfiltered.narrowPeak"), paramset=config['macs2_narrow_params'].keys(), sample=config['sample_paths'].keys()),
            expand(os.path.join(MACS2_NARROW_DIR, "{paramset}", "{sample}.summits.BLfiltered.bed"), paramset=config['macs2_narrow_params'].keys(), sample=config['sample_paths'].keys()),
        ]
    else:
        all_input = [
            expand(os.path.join(BIGWIG_DIR, "{sample}.1m.bw"), sample=config['sample_paths'].keys()), #Create bigwigs for all samples
            expand(os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.BLfiltered.hpeaks"), paramset=config['homer_findPeaks_params'].keys(), sample=config['sample_input'].keys()), #Call peaks for all samples with matched inputs
            expand(os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.BLfiltered.narrowPeak"), paramset=config['macs2_broad_params'].keys(), sample=config['sample_paths'].keys()),
            expand(os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.summits.BLfiltered.bed"), paramset=config['macs2_broad_params'].keys(), sample=config['sample_paths'].keys()),
            expand(os.path.join(MACS2_NARROW_DIR, "{paramset}", "{sample}.BLfiltered.narrowPeak"), paramset=config['macs2_narrow_params'].keys(), sample=config['sample_paths'].keys()),
            expand(os.path.join(MACS2_NARROW_DIR, "{paramset}", "{sample}.summits.BLfiltered.bed"), paramset=config['macs2_narrow_params'].keys(), sample=config['sample_paths'].keys()),
        ]
else:
    all_input = [
        expand(os.path.join(BIGWIG_DIR, "{sample}.1m.bw"), sample=config['sample_paths'].keys()), #Create bigwigs for all samples
        expand(os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.BLfiltered.narrowPeak"), paramset=config['macs2_broad_params'].keys(), sample=config['sample_paths'].keys()),
        expand(os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.summits.BLfiltered.bed"), paramset=config['macs2_broad_params'].keys(), sample=config['sample_paths'].keys()),
        expand(os.path.join(MACS2_NARROW_DIR, "{paramset}", "{sample}.BLfiltered.narrowPeak"), paramset=config['macs2_narrow_params'].keys(), sample=config['sample_paths'].keys()),
        expand(os.path.join(MACS2_NARROW_DIR, "{paramset}", "{sample}.summits.BLfiltered.bed"), paramset=config['macs2_narrow_params'].keys(), sample=config['sample_paths'].keys()),
    ]

# Print pipeline version number to log
include: "version.smk"
version_string = "Pipeline version: {}\n".format(version)
logger.logger.info(version_string)

# Rules

rule all:
    input:
        all_input

include: "rules/alignment_bwa_aln_pe.smk"
include: "rules/homer.smk"
include: "rules/homer_cutandrun.smk"

rule concatenate_reads:
    input:
        lambda wildcards: config['sample_paths'][wildcards.sample][wildcards.read]
    output:
        os.path.join(CONCAT_READS_DIR, "{sample}_R{read}.fastq.gz")
    shell:
        "cat {input} > {output}"

rule mark_duplicates:
    input:
        os.path.join(ALIGN_DIR, "{sample}.sorted.bam")
    output:
        bam = os.path.join(ALIGN_DIR, "{sample}.mrkdup.bam"),
        metric = os.path.join(ALIGN_DIR, "{sample}.mrkdup.metric"),
        condaenv = os.path.join(CONDA_ENVS_DIR, "{sample}.picard.env.txt")
    params:
        tmpdir = config['tmpdir']
    conda: "envs/picard.yaml"
    shell:
        "export JAVA_OPTIONS=-Xmx12g ; "
        "picard MarkDuplicates I={input} O={output.bam} "
        "METRICS_FILE={output.metric} "
        "ASSUME_SORTED=True "
        "VALIDATION_STRINGENCY=LENIENT "
        "TMP_DIR={params.tmpdir} ;"
        "conda list --export > {output.condaenv}"

# samtools is available in the parent environment atac_chip_pipeline
rule index_dupmarked_bams:
    input:
        os.path.join(ALIGN_DIR, "{sample}.mrkdup.bam")
    output:
        os.path.join(ALIGN_DIR, "{sample}.mrkdup.bai")
    shell:
        "samtools index {input} {output}"

# samtools is available in the parent environment atac_chip_pipeline
rule samtools_prune:
    input:
        bam = os.path.join(ALIGN_DIR, "{sample}.mrkdup.bam"),
        bai = os.path.join(ALIGN_DIR, "{sample}.mrkdup.bai")
    output:
        bam = temp(os.path.join(PRUNE_DIR, "{sample}.stpruned.bam"))
    params:
        incl_chr = lambda wildcards: INCLUDE_CHRS[get_genome(wildcards.sample)],
        flags = config['samtools_prune_flags']
    shell:
        "samtools view -b {params.flags} {input.bam} {params.incl_chr} > {output.bam}"

# samtools is available in the parent environment atac_chip_pipeline
rule namesort_st_pruned:
    input:
        os.path.join(PRUNE_DIR, "{sample}.stpruned.bam")
    output:
        temp(os.path.join(PRUNE_DIR, "{sample}.ns.bam"))
    shell:
        "samtools sort -n -o {output} {input}"

rule X0_pair_filter:
    input:
        os.path.join(PRUNE_DIR, "{sample}.ns.bam")
    output:
        filtered = temp(os.path.join(PRUNE_DIR, "{sample}.x0_filtered.bam")),
        condaenv = os.path.join(CONDA_ENVS_DIR, "{sample}.pysam.env.txt")
    params:
        config['X0_pair_filter_params']
    conda: "envs/pysam.yaml"
    shell:
        "python {SCRIPTS_DIR}/X0_pair_filter.py {params} -b {input} -o {output.filtered} ;"
        "conda list --export > {output.condaenv}"

# samtools is available in the parent environment atac_chip_pipeline
rule coordsort_index_final_pruned:
    input:
        os.path.join(PRUNE_DIR, "{sample}.x0_filtered.bam")
    output:
        bam = os.path.join(PRUNE_DIR, "{sample}.pruned.bam"),
        bai = os.path.join(PRUNE_DIR, "{sample}.pruned.bai")
    shell:
        "samtools sort -o {output.bam} {input} ;"
        "samtools index {output.bam} {output.bai}"

rule deeptools_bamcoverage_bw:
    input:
        bam = os.path.join(PRUNE_DIR, "{sample}.pruned.bam"),
        bai = os.path.join(PRUNE_DIR, "{sample}.pruned.bai")
    output:
        bigwig = os.path.join(BIGWIG_DIR, "{sample}.1m.bw"),
        condaenv = os.path.join(CONDA_ENVS_DIR, "{sample}.deeptools.env.txt")
    params:
        blacklist = lambda wildcards: "-bl {}".format(config['blacklist'][get_genome(wildcards.sample)]) if config.get('deeptools_bamcoverage_use_blacklist') == True else '', # deeptools_bamcoverage_use_blacklist should be True or False. If True, add `-bl /path/to/blacklist` to command. Otherwise empty string
        args = config['deeptools_bamcoverage_params']
    conda: "envs/deeptools.yaml"
    shell:
        "bamCoverage --bam {input.bam} -o {output.bigwig} {params.blacklist} {params.args}; "
        "conda list --export > {output.condaenv}"

rule deeptools_bamcoverage_bg:
    input:
        bam = os.path.join(PRUNE_DIR, "{sample}.pruned.bam"),
        bai = os.path.join(PRUNE_DIR, "{sample}.pruned.bai")
    output:
        os.path.join(BEDGRAPH_DIR, "{sample}.1m.bedgraph")
    params:
        blacklist = lambda wildcards: "-bl {}".format(config['blacklist'][get_genome(wildcards.sample)]) if config.get('deeptools_bamcoverage_use_blacklist') == True else '', # deeptools_bamcoverage_use_blacklist should be True or False. If True, add `-bl /path/to/blacklist` to command. Otherwise empty string
        args = config['deeptools_bamcoverage_params']
    conda: "envs/deeptools.yaml"
    shell:
        "bamCoverage --bam {input.bam} --outFileFormat bedgraph -o {output} {params.blacklist} {params.args}"

rule macs2_narrow_peaks:
    input:
        sample = os.path.join(PRUNE_DIR, "{sample}.pruned.bam"),
        input = lambda wildcards: os.path.join(HOMERTAG_DIR, config['sample_input'][wildcards.sample], ".pruned.bam")
    output:
        temp(os.path.join(MACS2_NARROW_DIR, "{paramset}", "{sample}.peaks.narrowPeak")),
        os.path.join(MACS2_NARROW_DIR, "{paramset}", "{sample}.peaks.xls"),
        temp(os.path.join(MACS2_NARROW_DIR, "{paramset}", "{sample}.summits.bed")),
        temp(os.path.join(MACS2_NARROW_DIR, "{paramset}", "{sample}.macs2.env.txt"))
    params:
        name = "{sample}",
        genome = lambda wildcards: MACS2_GENOME_SIZE[get_genome(wildcards.sample)],
        outdir = os.path.join(MACS2_NARROW_DIR, "{paramset}"),
        macs2_narrow_params = lambda wildcards: config['macs2_narrow_params'][wildcards.paramset],
        condadir = os.path.join(CONDA_ENVS_DIR)
    conda: "envs/macs2.yaml"
    run:
      if input.sample == input.input:
          shell("macs2 callpeak -t {input.sample} --outdir {params.outdir} -n {params.name} -g {params.genome} {params.macs2_narrow_params} ; ")
          # aligning filename formats - separated by '.'s preferred
          shell("mv {params.outdir}/{wildcards.sample}_peaks.xls {params.outdir}/{wildcards.sample}.peaks.xls ; ")
          shell("mv {params.outdir}/{wildcards.sample}_peaks.narrowPeak {params.outdir}/{wildcards.sample}.peaks.narrowPeak ; ")
          shell("mv {params.outdir}/{wildcards.sample}_summits.bed {params.outdir}/{wildcards.sample}.summits.bed ;")
          shell("conda list --export > {wildcards.sample}.macs2.env.txt ;")
          shell("mv {params.outdir}/{wildcards.sample}.macs2.env.txt {params.condadir}/{wildcards.sample}.macs2.env.txt ;")
      else:
          shell("macs2 callpeak -t {input.sample} -c {input.input} --outdir {params.outdir} -n {params.name} -g {params.genome} {params.macs2_narrow_params} ; ")
          # aligning filename formats - separated by '.'s preferred
          shell("mv {params.outdir}/{wildcards.sample}_peaks.xls {params.outdir}/{wildcards.sample}.peaks.xls ; ")
          shell("mv {params.outdir}/{wildcards.sample}_peaks.narrowPeak {params.outdir}/{wildcards.sample}.peaks.narrowPeak ; ")
          shell("mv {params.outdir}/{wildcards.sample}_summits.bed {params.outdir}/{wildcards.sample}.summits.bed ;")
          shell("conda list --export > {wildcards.sample}.macs2.env.txt ;")
          shell("mv {params.outdir}/{wildcards.sample}.macs2.env.txt {params.condadir}/{wildcards.sample}.macs2.env.txt ;")

rule macs2_broad_peaks:
    input:
        sample = os.path.join(PRUNE_DIR, "{sample}.pruned.bam"),
        input = lambda wildcards: os.path.join(HOMERTAG_DIR, config['sample_input'][wildcards.sample], ".pruned.bam")
    output:
        temp(os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.peaks.narrowPeak")),
        os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.peaks.xls"),
        temp(os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.summits.bed")),
        temp(os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.peaks.broadPeak")),
        temp(os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.peaks.gappedPeak"))
    params:
        name = "{sample}",
        genome = lambda wildcards: MACS2_GENOME_SIZE[get_genome(wildcards.sample)],
        outdir = os.path.join(MACS2_BROAD_DIR, "{paramset}"),
        macs2_broad_params = lambda wildcards: config['macs2_broad_params'][wildcards.paramset]
    conda: "envs/macs2.yaml"
    run:
      if input.sample == input.input:
          shell("macs2 callpeak -t {input.sample} --outdir {params.outdir} -n {params.name} -g {params.genome} {params.macs2_narrow_params} ; ")
          # aligning filename formats - separated by '.'s preferred
          shell("mv {params.outdir}/{wildcards.sample}_peaks.xls {params.outdir}/{wildcards.sample}.peaks.xls ; ")
          shell("mv {params.outdir}/{wildcards.sample}_peaks.narrowPeak {params.outdir}/{wildcards.sample}.peaks.narrowPeak ; ")
          shell("mv {params.outdir}/{wildcards.sample}_summits.bed {params.outdir}/{wildcards.sample}.summits.bed ;")
          shell("mv {params.outdir}/{wildcards.sample}_peaks.broadPeak {params.outdir}/{wildcards.sample}.peaks.broadPeak ; ")
          shell("mv {params.outdir}/{wildcards.sample}_peaks.gappedPeak {params.outdir}/{wildcards.sample}.peaks.gappedPeak ; ")
      else:
          shell("macs2 callpeak -t {input.sample} -c {input.input} --outdir {params.outdir} -n {params.name} -g {params.genome} {params.macs2_narrow_params} ; ")
          # aligning filename formats - separated by '.'s preferred
          shell("mv {params.outdir}/{wildcards.sample}_peaks.xls {params.outdir}/{wildcards.sample}.peaks.xls ; ")
          shell("mv {params.outdir}/{wildcards.sample}_peaks.narrowPeak {params.outdir}/{wildcards.sample}.peaks.narrowPeak ; ")
          shell("mv {params.outdir}/{wildcards.sample}_summits.bed {params.outdir}/{wildcards.sample}.summits.bed ;")
          shell("mv {params.outdir}/{wildcards.sample}_peaks.broadPeak {params.outdir}/{wildcards.sample}.peaks.broadPeak ; ")
          shell("mv {params.outdir}/{wildcards.sample}_peaks.gappedPeak {params.outdir}/{wildcards.sample}.peaks.gappedPeak ; ")

rule blacklist_filter_narrow_peaks:
    input:
        narrowpeak = os.path.join(MACS2_NARROW_DIR, "{paramset}", "{sample}.peaks.narrowPeak"),
        summits = os.path.join(MACS2_NARROW_DIR, "{paramset}", "{sample}.summits.bed")
    output:
        narrowpeak = os.path.join(MACS2_NARROW_DIR, "{paramset}", "{sample}.BLfiltered.narrowPeak"),
        summits = os.path.join(MACS2_NARROW_DIR, "{paramset}", "{sample}.summits.BLfiltered.bed"),
        condaenv = os.path.join(CONDA_ENVS_DIR, "{sample}.bedtools.env.txt")
    params:
        blacklist = lambda wildcards: config['blacklist'][get_genome(wildcards.sample)]
    conda: "envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.narrowpeak} -b {params.blacklist} -v > {output.narrowpeak} ; "
        "bedtools intersect -a {input.summits} -b {params.blacklist} -v > {output.summits} ; "
        "conda list --export > {output.condaenv}"

rule blacklist_filter_broad_peaks:
    input:
        narrowpeak = os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.peaks.narrowPeak"),
        broadpeak = os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.peaks.broadPeak"),
        gappedpeak = os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.peaks.gappedPeak"),
        summits = os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.summits.bed")
    output:
        narrowpeak = os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.BLfiltered.narrowPeak"),
        broadpeak = os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.BLfiltered.broadPeak"),
        gappedpeak = os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.BLfiltered.gappedPeak"),
        summits = os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.summits.BLfiltered.bed")
    params:
        blacklist = lambda wildcards: config['blacklist'][get_genome(wildcards.sample)]
    conda: "envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.narrowpeak} -b {params.blacklist} -v > {output.narrowpeak} ; "
        "bedtools intersect -a {input.broadpeak} -b {params.blacklist} -v > {output.broadpeak} ; "
        "bedtools intersect -a {input.gappedpeak} -b {params.blacklist} -v > {output.gappedpeak} ; "
        "bedtools intersect -a {input.summits} -b {params.blacklist} -v > {output.summits} ; "