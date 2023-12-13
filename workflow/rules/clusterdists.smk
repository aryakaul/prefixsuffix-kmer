checkpoint cluster_dists:
    output:
        directory(f"{dir_output()}" + "/passing_genes/{batch}/{genebatch}"),
    input:
        filterdists=aggregate_filter_dists,
    conda:
        "../envs/sklearn.yml"
    params:
        intermediate=config['intermediate_dir'],
        script=snakemake.workflow.srcdir("../scripts/cluster_dists"),
    shell:
        """
        mkdir -p {output}
        {params.script} \\
                -i {params.intermediate}/kmerdists/{wildcards.batch}/{wildcards.genebatch} \\
                -o {output} \\
                -vv \\
                --epsilon 500 \\
                --min_samples 1 \\

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
