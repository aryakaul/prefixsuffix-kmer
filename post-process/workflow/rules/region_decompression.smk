checkpoint region_decompression:
    input:
        downsamplecsv=fn_downsampled_df("{batch}", "{genebatch}", "{passinggene}"),
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
        kmerlen=config["kmer_length"],
        script=Path(workflow.basedir) / "scripts/region_decompression",
    shell:
        """
        {params.script} \\
            -i {input.downsamplecsv} \\
            -f {input.fof} \\
            -e {output} \\
            -w {params.windowsize} \\
            -m \\
            -k {params.kmerlen} \\
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
        
        # First, escape dots. Note the double escaping: once for the shell and once for sed.
        # Second, escape plus signs.
        # Third, replace ':' with '.', which doesn't need escaping in this context.
        regex=$(echo "$x" | sed -e 's/\\./\\\\./g' -e 's/\\+/\\\\+/g' -e 's/:/./g')
        
        seqkit grep -n -r -p "$regex" {input} > {output}
        """


rule aggregate_regionalfastas:
    output:
        fn_decompregion_agg("{batch}", "{genebatch}", "{passinggene}"),
    input:
        #chkpntaggregate_regionalfastas,
        chkpntaggregate_blastouts,
        #fn_blastouts("{batch}", "{genebatch}", "{passinggene}"),
        #expand(
            #fn_blastrawout("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
            #contig=glob_wildcards("{input.reffastadir}/{contig}.fasta").contig
        #)
    shell:
        """
        echo {input} > {output}
        """
