rule bakta_annotation:
    input:
        fa=f"{dir_prevoutput()}/decompressed_genomes/genomes" + "/{batch}/{genome}.fa",
    output:
        genbank=fn_bakta_gff("{batch}", "{genome}"),
        faa=fn_bakta_faa("{batch}", "{genome}"),
    conda:
        "../envs/bakta.yml"
    params:
        output_dir=f"{dir_intermediate()}/bakta_out" + "/{batch}/{genome}",
        prefix="{genome}",
        bakta_db=config["bakta_db"],
    shell:
        """
        bakta \\
            --db {params.bakta_db}  \\
            --output {params.output_dir}  \\
            --prefix {params.prefix}  \\
            --force \\
            --locus-tag {params.prefix} \\
            {input.fa}
        """
