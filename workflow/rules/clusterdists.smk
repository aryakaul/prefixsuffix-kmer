checkpoint cluster_dists:
    output:
        directory(f"{dir_output()}" + "/passing_genes/{batch}/{genebatch}"),
    input:
        filterdists=aggregate_filter_dists,
    conda:
        "../envs/sklearn.yml"
    params:
        intermediate=config['intermediate_dir'],
        min_samples=config['min_samples'],
        min_cluster_size=config['min_cluster_size'],
        cluster_selection_epsilon=config['cluster_selection_epsilon'],
        metric=config['metric'],
        script=snakemake.workflow.srcdir("../scripts/cluster_dists"),
    shell:
        """
        mkdir -p {output}
        {params.script} \\
                -i {params.intermediate}/kmerdists/{wildcards.batch}/{wildcards.genebatch} \\
                -o {output} \\
                -vv \\
                --min_cluster_size {params.min_cluster_size} \\
                --min_samples {params.min_samples} \\
                --cluster_selection_epsilon {params.cluster_selection_epsilon} \\
                --metric {params.metric}
        """

rule aggregate_2:
    input:
        aggregate_passing_genes,
    output:
        f"{dir_output()}" + "/passing_genes-{batch}-{genebatch}.txt"
    shell:
        """
        echo {input} > {output}
        """
