rule fastmap_process:
    output:
        fn_fastmapprocess(
            _batch="{batch}", _bucket="{bucket}", _genebatch="{genebatch}"
        ),
    input:
        fn_fastmapraw(_batch="{batch}", _bucket="{bucket}", _genebatch="{genebatch}"),
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


rule fastmap_distances:
    output:
        fn_fastmapdists(_batch="{batch}", _bucket="{bucket}", _genebatch="{genebatch}"),
    input:
        fn_fastmapprocess(
            _batch="{batch}", _bucket="{bucket}", _genebatch="{genebatch}"
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


rule aggregate:
    input:
        aggregate_fastmap_dists,
    output:
        f"{dir_output()}" + "/{batch}-fastmapdists.txt",
    shell:
        """
        echo {input} > {output}
        """
