# def get_time(wildcards, attempt):
# if attempt == 1:
# return "0-02:00:00"
# if attempt == 2:
# return "0-05:00:00"
# if attempt == 3:
# return "0-10:00:00"


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
        script=Path(workflow.basedir) / "scripts/process_fastmapout",
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
        script=Path(workflow.basedir) / "scripts/analyze_fastmapout",
        kmer_length=config["kmer_length"],
        gap_distance=config["gap_distance"],
    resources:
        mem_gb=10,
    shell:
        """
        {params.script} \\
            -i {input} \\
            -o {output} \\
            -k {params.kmer_length} \\
            -g {params.gap_distance} \\
            -vv
        """


rule parse_distances:
    output:
        fn_parsedists(_batch="{batch}", _bucket="{bucket}", _genebatch="{genebatch}"),
    input:
        fn_fastmapdists(_batch="{batch}", _bucket="{bucket}", _genebatch="{genebatch}"),
    conda:
        "../envs/pandas.yml"
    params:
        script=Path(workflow.basedir) / "scripts/parse_distances",
    # resources:
    # time=get_time
    shell:
        """
        mkdir -p $(dirname {output})
        {params.script} \\
        {input} \\
        {output} 
        """


rule aggregate:
    input:
        aggregate_filter_dists,
    output:
        f"{dir_output()}" + "/{batch}-filterdists.txt",
    shell:
        """
        echo {input} > {output}
        """
