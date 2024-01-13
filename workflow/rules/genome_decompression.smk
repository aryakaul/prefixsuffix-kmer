rule passinggene_cluster_decompression:
    output:
        clustersampled_csv=fn_downsampled_df("{batch}", "{genebatch}", "{passing_gene}"),
    input:
        #passing_gene_csv=f"{dir_output()}" + "/passing_genes-{batch}-{genebatch}.txt",
        passing_gene_csv=f"{dir_output()}"
        + "/{batch}/{genebatch}/passing_genes/{passing_gene}/{passing_gene}_clusters.csv",
        fof=f"{dir_input()}" + "/{batch}.txt",
    params:
        script=Path(workflow.basedir) / "scripts/genome_decompression",
        sample_size=1,
    shell:
        """
        OUTDIR=$(dirname {output})/genomes
        {params.script} \\
            -i {input.passing_gene_csv} \\
            -f {input.fof} \\
            -e $OUTDIR \\
            -o {output.clustersampled_csv} \\
            -s {params.sample_size} \\
            -vvv 
        """
