checkpoint region_decompression:
    input:
        clustercsv=fn_cluster_csv("{batch}", "{genebatch}", "{passinggene}"),
        fof=fn_fof("{batch}"),
    output:
        directory(
            f"{dir_intermediate()}"
            + "/decompressed_regions/{batch}/{genebatch}/{passinggene}"
        ),
    conda:
        "../envs/biopython.yml"
    params:
        windowsize=config["window_size"],
        script=Path(workflow.basedir) / "scripts/region_decompression",
    shell:
        """
        {params.script} \\
            -i {input.clustercsv} \\
            -f {input.fof} \\
            -e {output} \\
            -w {params.windowsize} \\
            -vvv
        """


rule passinggene_fasta:
    input:
        fn_genebatch_input("{genebatch}"),
    output:
        fn_passinggenefasta("{batch}", "{genebatch}", "{passinggene}"),
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        x=$(echo "{wildcards.passinggene}" | sed 's/-len.[0-9]*//')
        seqkit grep -n -r -p "$x" {input} > {output}
        """


rule region_minimap:
    input:
        reffasta=fn_regionfa("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
        queryfasta=fn_passinggenefasta("{batch}", "{genebatch}", "{passinggene}"),
    output:
        fn_minimaprawout("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
    conda:
        "../envs/minimap2.yml"
    shell:
        """
        minimap2 -x splice -un -k8 {input.reffasta} {input.queryfasta} > {output}
        """


rule aggregate_regionalfastas:
    output:
        fn_decompregion_agg("{batch}", "{genebatch}", "{passinggene}"),
    input:
        chkpntaggregate_regionalfastas,
    shell:
        """
        echo {input} > {output}
        """
