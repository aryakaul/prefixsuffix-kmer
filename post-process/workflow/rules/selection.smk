rule selection_analysis:
    output:
        fn_mutationjson("{batch}", "{genebatch}", "{passinggene}"),
        fn_dndsanalysis("{batch}", "{genebatch}", "{passinggene}"),
    input:
        fastadir=chkpntaggregate_regionalfasta_dir,
        knownorf=fn_passinggenefasta("{batch}", "{genebatch}", "{passinggene}"),
    params:
        min_orf_len=100,
        cdhit_identity=0.9,
        cdhit_length=0.8,
        script=Path(workflow.basedir) / "scripts/selection_analysis",
    threads: 16
    conda:
        "../envs/pyrodigal.yml"
    shell:
        """
        OUTPUT=$(dirname {output.fn_mutationjson})
        {params.script} \\
            {input.fastadir} \\
            \$OUTPUT \\
            --min_orf_len {params.min_orf_len} \\
            --cdhit_identity {params.cdhit_identity} \\
            --cdhit_length {params.cdhit_length} \\
            --threads {threads} \\
            --known_orf_fasta {input.knownorf} \\
            -vv
        """
