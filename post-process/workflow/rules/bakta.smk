rule bakta_annotation:
    input:
        downsampled_annot_df=fn_tinydownsampled_df(
            "{batch}", "{genebatch}", "{passinggene}"
        ),
    output:
        bakta_annot_flag=fn_bakta_annot_done("{batch}", "{genebatch}", "{passinggene}"),
    conda:
        "../envs/bakta.yml"
    params:
        script=Path(workflow.basedir) / "scripts/download_annotation",
    shell:
        """
        mkdir -p $(dirname {output})
        {params.script} {input} $(dirname {output})
        touch {output}
        """
