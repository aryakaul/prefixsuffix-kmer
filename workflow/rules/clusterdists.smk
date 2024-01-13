checkpoint cluster_dists:
    output:
        directory(f"{dir_output()}" + "/{batch}/{genebatch}"),
    input:
        filterdists=aggregate_filter_dists,
    conda:
        "../envs/sklearn.yml"
    params:
        epsilon=config["max_distance"],
        min_samples=config["min_samples"],
        intermediate=config["intermediate_dir"],
        script=Path(workflow.basedir) / "scripts/cluster_dists",
    shell:
        """
        mkdir -p {output}
        {params.script} \\
                -i {params.intermediate}/kmerdists/{wildcards.batch}/{wildcards.genebatch} \\
                -o {output} \\
                --epsilon {params.epsilon} \\
                --min_samples {params.min_samples} \\
                -vv 
        """


rule aggregate_2:
    output:
        f"{dir_output()}" + "/passing_genes-{batch}-{genebatch}.txt",
    input:
        aggregate_passing_genes,
    shell:
        """
        echo {input} > {output}
        """
