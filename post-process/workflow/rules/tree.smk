rule build_sketched_tree:
    output:
        tree_newick=fn_downsample_nwk("{batch}", "{genebatch}", "{passinggene}"),
    input:
        chkpntaggregate_genomefasta_dir,
        filelist=fn_decompgenome_list("{batch}", "{genebatch}", "{passinggene}"),
    conda:
        "../envs/attotree.yml"
    shell:
        """
        attotree \\
            -L {input.filelist} \\
            -o {output.tree_newick} \\
            -D
        """


rule itol_annottext:
    output:
        colorrange=fn_colorrange_annot("{batch}", "{genebatch}", "{passinggene}"),
        simplebar=fn_simplebar_annot("{batch}", "{genebatch}", "{passinggene}"),
        newtree=fn_itoltree("{batch}", "{genebatch}", "{passinggene}"),
    input:
        tree_newick=fn_downsample_nwk("{batch}", "{genebatch}", "{passinggene}"),
        clustercsv=fn_downsampled_df("{batch}", "{genebatch}", "{passinggene}"),
    conda:
        "../envs/pandas.yml"
    params:
        metadata_661k=f"{dir_intermediate()}"
        + "/../supp_data/File2_taxid_lineage_661K.txt",
        script=Path(workflow.basedir) / "scripts/make_itol_annots",
    shell:
        """
        #if [ "{wildcards.batch}" == "661k_allsamples" ]; then
        #    {params.script} \\
                #        -i {input.clustercsv} \\
                #-t {input.tree_newick} \\
                #-m {params.metadata_661k} \\
                #-o $(dirname {output.newtree})
        #else
            {params.script} \\
                -i {input.clustercsv} \\
                -t {input.tree_newick} \\
                -o $(dirname {output.newtree})
        """


rule minimal_cuts:
    output:
        fn_mincuts("{batch}", "{genebatch}", "{passinggene}"),
    input:
        tree=fn_itoltree("{batch}", "{genebatch}", "{passinggene}"),
        csv=fn_downsampled_df("{batch}", "{genebatch}", "{passinggene}"),
    conda:
        "../envs/ete4.yml"
    params:
        script=Path(workflow.basedir) / "scripts/minimal_cuts",
    shell:
        """
        clusterzero=$(dirname {output})/$(basename {output})-cluster0
        tree=$(dirname {output})/$(basename {output})-tree
        csvtk filter -f 'Cluster=0' {input.csv} | csvtk cut -f GenomeID \\
            | sed 1d | csvtk transpose | csvtk csv2tab > $clusterzero
        {params.script} \\
            {input.tree} \\
            $clusterzero \\
            -t $tree \\
            -o {output} -v
        rm $tree
        """
