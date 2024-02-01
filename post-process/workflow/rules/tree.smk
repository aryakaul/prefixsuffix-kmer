rule build_sketched_tree:
    output:
        tree_newick=fn_downsample_nwk("{batch}", "{genebatch}", "{passinggene}"),
        filelist=f"{dir_intermediate()}"
        + "/decompressed_genomes/{batch}/{genebatch}/{passinggene}_filelist.txt",
    input:
        chkpntaggregate_genomefasta_dir,
    params:
        script=Path(workflow.basedir) / "scripts/attotree.py",
    conda:
        "../envs/attotree.yml"
    shell:
        """
        for i in {input}/*.fa; do echo $i; done > {output.filelist} 

        {params.script} \\
            -L {output.filelist} \\
            -o {output.tree_newick}
        """


rule itol_annottext:
    output:
        colorrange=fn_colorrange_annot("{batch}", "{genebatch}", "{passinggene}"),
        simplebar=fn_simplebar_annot("{batch}", "{genebatch}", "{passinggene}"),
        newtree=fn_newtree("{batch}", "{genebatch}", "{passinggene}"),
    input:
        tree_newick=fn_downsample_nwk("{batch}", "{genebatch}", "{passinggene}"),
    conda:
        "../envs/pandas.yml"
    params:
        downsampled_csv=f"{dir_intermediate()}"
        + "/decompressed_genomes/{batch}/{genebatch}/{passinggene}_downsampled.csv",
        script=Path(workflow.basedir) / "scripts/make_itol_annots",
    shell:
        """
        {params.script} \\
            -i {params.downsampled_csv} \\
            -t {input.tree_newick} \\
            -o $(dirname {output.newtree})
        """
