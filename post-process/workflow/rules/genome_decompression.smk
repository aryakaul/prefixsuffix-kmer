checkpoint genome_decompression:
    input:
        clustercsv=fn_downsampled_df("{batch}", "{genebatch}", "{passinggene}"),
        fof=fn_fof("{batch}"),
    output:
        outdir=directory(
            f"{dir_intermediate()}"
            + "/decompressed_genomes/{batch}/{genebatch}/{passinggene}"
        ),
    conda:
        "../envs/pandas.yml"
    params:
        script=Path(workflow.basedir) / "scripts/genome_decompression",
    shell:
        """
        {params.script} \\
            -i {input.clustercsv} \\
            -f {input.fof} \\
            -e {output.outdir} \\
            -vvv
        """


rule aggregate_fullgenomefasta:
    output:
        fn_decompgenome_agg("{batch}", "{genebatch}", "{passinggene}"),
    input:
        chkpntaggregate_genomefastas,
    shell:
        """
        echo {input} > {output}
        """
