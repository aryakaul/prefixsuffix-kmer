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


rule make_blastidx:
    input:
        reffasta=fn_regionfa("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
    output:
        fn_blastidx("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
    conda:
        "../envs/blast.yml"
    shell:
        """
        makeblastdb -in {input} -dbtype nucl
        """


checkpoint run_blast:
    input:
        reffastas=chkpntaggregate_regionalfastas,
        #reffastadir=fn_regionfadir("{batch}", "{genebatch}", "{passinggene}"),
        #reffasta=fn_regionfa("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
        queryfasta=fn_passinggenefasta("{batch}", "{genebatch}", "{passinggene}"),
    output:
        #chkpntaggregate_blastouts,
        #fn_blastouts("{batch}", "{genebatch}", "{passinggene}")
        #expand(
        #fn_blastrawout("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
        #contig=glob_wildcards("{input.reffastadir}/{contig}.fasta").contig
        #)
        #expand(
        directory(
            f"{dir_intermediate()}" + "/blast/raw/{batch}/{genebatch}/{passinggene}"
        ),
        #fn_blastrawout("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
        #contig=wildcards.contig,
        #)
    conda:
        "../envs/blast.yml"
    shell:
        #"""
        #blastn -query {input.queryfasta} -subject {input.reffasta} -outfmt 6 > {output}
        #"""
        """
        x="{input.reffastas}"
        mkdir -p {output}
        for i in $x; do
            blastn -query {input.queryfasta} -subject $i -outfmt 6 \\
                    > {output}/$(basename $i .fasta).blastout
        done
        """
