checkpoint make_bwaidx:
    output:
        directory(f"{dir_intermediate()}" + "/bwa_indices/{batch}"),
    input:
        fof=f"{dir_input()}" + "/{batch}.txt",
    params:
        script=snakemake.workflow.srcdir("../scripts/build_bwaindex_fof"),
        samplesperidx=config["batch_size"],
    conda:
        "../envs/bwamem.yml"
    shell:
        """
        {params.script}  \\
                {input} \\
                {params.samplesperidx} \\
                -o {output}
        """


rule fastmap:
    output:
        fn_fastmapraw(_batch="{batch}", _bucket="{bucket}", _genebatch="{genebatch}"),
    input:
        prefsuffkmers=fn_prefsuffkmer(_genebatch="{genebatch}"),
        bwaindex=fn_bwaidxbucket(_batch="{batch}", _bucket="{bucket}"),
    params:
        kmer_length=config["kmer_length"],
    conda:
        "../envs/bwamem.yml"
    shell:
        """
        mkdir -p $(dirname {output})
        idx=$(dirname {input.bwaindex})/$(basename {input.bwaindex} .bwt)
        bwa fastmap -w 99999 -l {params.kmer_length} $idx {input.prefsuffkmers} > {output}
        """


# rule aggregate:
#    input:
#        aggregate_fastmap_raw,
#    output:
#        f"{dir_output()}" + "/{batch}-fastmapraw.txt",
#    shell:
#        "cat {input} > {output}"
