rule process_fastmap_out:
    output:
        fastmapprocessed=fn_fastmapprocess(
            _sample="{sample}", _batch="{batch}", _genebatch="{genebatch}"
        ),
    input:
        fastmapout=fn_fastmapraw(
            _sample="{sample}", _batch="{batch}", _genebatch="{genebatch}"
        ),
    conda:
        "../envs/biopython.yml"
    params:
        script=snakemake.workflow.srcdir("../scripts/process_fastmapout"),
    shell:
        """
        {params.script} \\
                -i {input} \\
                -o {output}
        """


rule analyze_fastmap_out:
    output:
        fastmapdistances=fn_fastmapdists(
            _sample="{sample}", _batch="{batch}", _genebatch="{genebatch}"
        ),
    input:
        fastmapprocessed=fn_fastmapprocess(
            _sample="{sample}", _batch="{batch}", _genebatch="{genebatch}"
        ),
    conda:
        "../envs/pandas.yml"
    params:
        script=snakemake.workflow.srcdir("../scripts/analyze_fastmapout"),
        kmer_length=config["kmer_length"],
        gap_distance=config["gap_distance"],
    shell:
        """
        {params.script} \\
                -i {input} \\
                -o {output} \\
                -k {params.kmer_length} \\
                -g {params.gap_distance}
        """
