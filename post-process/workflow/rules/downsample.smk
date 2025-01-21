rule downsample_df:
    input:
        clustercsv=fn_cluster_csv("{batch}", "{genebatch}", "{passinggene}"),
        fof=fn_fof("{batch}"),
    output:
        output_csv=fn_downsampled_df("{batch}", "{genebatch}", "{passinggene}"),
    conda:
        "../envs/pandas.yml"
    params:
        numbergenomes=config["number_genomes"],
        metadata661k=f"{dir_intermediate()}"
        + "/../supp_data/File4_QC_characterisation_661K.txt",
        metadataallthebact=f"{dir_intermediate()}"
        + "/../supp_data/hq_set.sample_list.txt.gz",
        script=Path(workflow.basedir) / "scripts/sample_genomes",
    shell:
        """
        mkdir -p $(dirname {output.output_csv})
        if [ "{wildcards.batch}" = "661k_allsamples" ]; then
            {params.script} \\
                -m {params.numbergenomes} \\
                -i {input.clustercsv} \\
                -o {output.output_csv} \\
                -q {params.metadata661k} \\
                -vvv
        elif [ "{wildcards.batch}" = "allthebacteriav2" ]; then
            {params.script} \\
                -m {params.numbergenomes} \\
                -i {input.clustercsv} \\
                -o {output.output_csv} \\
                -q {params.metadataallthebact} \\
                -a \\
                -vvv
        else
            {params.script} \\
                -m {params.numbergenomes} \\
                -i {input.clustercsv} \\
                -o {output.output_csv} \\
                -vvv
        fi
        """
