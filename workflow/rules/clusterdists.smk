checkpoint cluster_dists:
    output:
        directory(f"{dir_output()}" + "/{batch}/{genebatch}"),
    input:
        filterdists=aggregate_filter_dists,
    conda:
        "../envs/dbscan.yml"
    params:
        epsilon=config["max_distance"],
        min_samples=config["min_samples"],
        intermediate=config["intermediate_dir"],
        remove_outliers=config["remove_outliers"],
        script=Path(workflow.basedir) / "scripts/cluster_dists_duck",
    shell:
        """
        mkdir -p {output}
        {params.script} \\
                -i {params.intermediate}/kmerdists/{wildcards.batch}/{wildcards.genebatch} \\
                -o {output} \\
                --epsilon {params.epsilon} \\
                --min_samples {params.min_samples} \\
                --remove_outliers \\
                -vv 
        """


rule aggregate_2:
    output:
        f"{dir_output()}" + "/passing_genes-{batch}-{genebatch}.txt",
    input:
        aggregate_passing_genes,
    resources:
        mem_gb=1,
        time="0-00:05:00",
        threads=2,
    shell:
        """
        echo {input} > {output}
        """
