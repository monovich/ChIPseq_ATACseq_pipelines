rule findMotifsGenome_macs2_broad:
    input:
        os.path.join(MACS2_BROAD_DIR, "{paramset}", "{sample}.BLfiltered.narrowPeak")
    output:
        directory(os.path.join(MACS2_BROAD_HOMERMOTIF_DIR, "{paramset}", "{sample}"))
    params:
        genome = lambda wildcards: config['sample_homer_fmg_genome'][wildcards.sample],
        params = config['homer_fmg_params']
    shell:
        "findMotifsGenome.pl {input} {params.genome} {output} {params.params}"

rule findMotifsGenome_macs2_narrow:
    input:
        os.path.join(MACS2_NARROW_DIR, "{paramset}", "{sample}.BLfiltered.narrowPeak")
    output:
        directory(os.path.join(MACS2_NARROW_HOMERMOTIF_DIR, "{paramset}", "{sample}"))
    params:
        genome = lambda wildcards: config['sample_homer_fmg_genome'][wildcards.sample],
        params = config['homer_fmg_params']
    shell:
        "findMotifsGenome.pl {input} {params.genome} {output} {params.params}"
