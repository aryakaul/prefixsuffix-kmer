rule passinggene_cluster_decompression:
    output:
        clustersampled_csv=fn_downsampled_df("{batch}", "{genebatch}", "{passing_gene}"),
    input:
        passing_gene_csv=f"{dir_output()}"
        + "/{batch}/{genebatch}/passing_genes/{passing_gene}/{passing_gene}_clusters.csv",
        fof=f"{dir_input()}" + "/{batch}.txt",
    conda:
        "../envs/pandas.yml"
    params:
        script=Path(workflow.basedir) / "scripts/genome_decompression",
        sample_size=config['cluster_sample_size'],
        outdir=f"{dir_intermediate()}/decompressed_genomes/genomes" 
    shell:
        """
        {params.script} \\
            -i {input.passing_gene_csv} \\
            -f {input.fof} \\
            -e {params.outdir} \\
            -o {output.clustersampled_csv} \\
            -s {params.sample_size} \\
            -vvv 
        """
