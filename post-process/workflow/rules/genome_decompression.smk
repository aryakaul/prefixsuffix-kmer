checkpoint genome_decompression:
    input:
        clustercsv=fn_cluster_csv("{batch}", "{genebatch}", "{passinggene}"),
        fof=fn_fof("{batch}"),
    output:
        directory(
            f"{dir_intermediate()}"
            + "/decompressed_genomes/{batch}/{genebatch}/{passinggene}"
        ),
    conda:
        "../envs/pandas.yml"
    params:
        numbergenomes=config["number_genomes"],
        output_csv=f"{dir_intermediate()}"
        + "/decompressed_genomes/{batch}/{genebatch}/{passinggene}_downsampled.csv",
        metadata=f"{dir_intermediate()}"
        + "/../supp_data/File4_QC_characterisation_661K.txt",
        script=Path(workflow.basedir) / "scripts/genome_decompression",
    shell:
        """
        if [ "{wildcards.batch}" = "661k_allsamples" ]; then
            {params.script} \\
                -m {params.numbergenomes} \\
                -i {input.clustercsv} \\
                -f {input.fof} \\
                -e {output} \\
                -o {params.output_csv} \\
                -q {params.metadata} \\
                -vvv
        else
            {params.script} \\
                -m {params.numbergenomes} \\
                -i {input.clustercsv} \\
                -f {input.fof} \\
                -e {output} \\
                -o {params.output_csv} \\
                -vvv
        fi
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
