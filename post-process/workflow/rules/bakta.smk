rule bakta_annotation:
    input:
        filelist=fn_decompgenome_list("{batch}", "{genebatch}", "{passinggene}"),
        downsampled_annot_df=fn_downsampled_annot_df(
            "{batch}", "{genebatch}", "{passinggene}"
        ),
    output:
        bakta_annot_flag=fn_bakta_annot_done("{batch}", "{genebatch}", "{passinggene}"),
    conda:
        "../envs/bakta.yml"
    params:
        #output_dir=f"{dir_intermediate()}/bakta_out" + "/{batch}/{genome}",
        #prefix="{genome}",
        bakta_db=config["bakta_db"],
    shell:
        """
        mkdir -p $(dirname {output})
        zcat {input.downsampled_annot_df} | sed 1d | while read -r p; do
            seqname=$(echo $p | cut -d, -f2)
            inputfa=$(grep "$seqname" {input.filelist})
            bakta \\
                --db {params.bakta_db}  \\
                --output $(dirname {output})/$seqname  \\
                --prefix $seqname \\
                --force \\
                --locus-tag $seqname \\
                $inputfa 
        done
        touch {output}
        """
