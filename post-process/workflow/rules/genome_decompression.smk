checkpoint genome_decompression:
    input:
        clustercsv=fn_downsampled_df("{batch}", "{genebatch}", "{passinggene}"),
        fof=fn_fof("{batch}"),
    output:
        outdir=directory(
            f"{dir_intermediate()}"
            + "/decompressed/genomes/large_sample/{batch}/{genebatch}/{passinggene}"
        ),
    conda:
        "../envs/pandas.yml"
    params:
        script=Path(workflow.basedir) / "scripts/genome_decompression",
    shell:
        """
        {params.script} \\
            -i {input.clustercsv} \\
            -f {input.fof} \\
            -e {output.outdir} \\
            -vvv
        """

checkpoint tiny_genome_decompression:
    input:
        clustercsv=fn_tinydownsampled_df("{batch}", "{genebatch}", "{passinggene}"),
        fof=fn_fof("{batch}"),
    output:
        outdir=directory(
            f"{dir_intermediate()}"
            + "/decompressed/genomes/tiny_sample/{batch}/{genebatch}/{passinggene}"
        ),
    conda:
        "../envs/pandas.yml"
    params:
        script=Path(workflow.basedir) / "scripts/genome_decompression",
    shell:
        """
        {params.script} \\
            -i {input.clustercsv} \\
            -f {input.fof} \\
            -e {output.outdir} \\
            -vvv
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


rule build_filelist:
    output:
        filelist=fn_decompgenome_list("{batch}", "{genebatch}", "{passinggene}"),
    input:
        chkpntaggregate_genomefasta_dir,
    shell:
        r"""
        set -euo pipefail
        shopt -s nullglob

        fasta_files=( {input}/*.fa )
        # If none, that's an error
        if [ ${{#fasta_files[@]}} -eq 0 ]; then
            echo "ERROR: no .fa files found in {input}" >&2
            exit 1
        fi

        # Otherwise write the list
        printf "%s\n" "${{fasta_files[@]}}" > {output.filelist}
        """

rule build_tinyfilelist:
    output:
        filelist=fn_tinydecompgenome_list("{batch}", "{genebatch}", "{passinggene}"),
    input:
        chkpntaggregate_tinygenomefasta_dir,
    shell:
        r"""
        set -euo pipefail
        shopt -s nullglob

        fasta_files=( {input}/*.fa )
        # If none, that's an error
        if [ ${{#fasta_files[@]}} -eq 0 ]; then
            echo "ERROR: no .fa files found in {input}" >&2
            exit 1
        fi

        # Otherwise write the list
        printf "%s\n" "${{fasta_files[@]}}" > {output.filelist}
        """
