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
        metadataallthebact=Path(workflow.basedir) / "../supp_data/hq_set.sample_list.txt.gz",
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

rule downsample_annotdf:
    input:
        prevdownsampled_df=fn_downsampled_df("{batch}", "{genebatch}", "{passinggene}"),
    output:
        downsampled_annot_df=fn_downsampled_annot_df("{batch}", "{genebatch}", "{passinggene}"),
    conda:
        "../envs/ete3.yml"
    params:
        numbergenomes=config["number_annotate"],
        metadata661k=f"{dir_intermediate()}"
        + "/../supp_data/File4_QC_characterisation_661K.txt",
        metadataallthebact=Path(workflow.basedir) / "../supp_data/hq_set.sample_list.txt.gz",
        script=Path(workflow.basedir) / "scripts/sample_genomes",
    shell:
        """
        # Create a temporary directory for the split files
        rm -rf $(dirname {output})/tmp_clusters

        mkdir -p $(dirname {output})/tmp_clusters

        # Split the file by 'Cluster'
        zcat {input} | csvtk split -f Cluster -o $(dirname {output})/tmp_clusters

        # For each split file, sample X random lines
        for f in $(dirname {output})/tmp_clusters/*.csv; do 
            {{  head -n 1 $f; tail -n +2 $f | shuf -n {params.numbergenomes}; }} > $(dirname $f)/$(basename $f .csv).csv_sample
        done

        # Concatenate all the sampled files into one CSV
        csvtk concat $(dirname {output})/tmp_clusters/*.csv_sample | gzip - > {output}
        rm -rf $(dirname {output})/tmp_clusters
        """
