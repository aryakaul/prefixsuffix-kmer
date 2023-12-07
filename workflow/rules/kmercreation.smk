rule prefixsuffix_kmergen:
    output:
        prefsuffkmerfasta=fn_prefsuffkmer(_genebatch="{genebatch}"),
    input:
        genebatch_ffn=f"{dir_input()}/{{genebatch}}.ffn",
    conda:
        "../envs/biopython.yml"
    params:
        script=snakemake.workflow.srcdir("../scripts/prefsuff_kmerextraction"),
        kmer_length=config["kmer_length"],
        gap_dist=config["gap_distance"],
    shell:
        """
        {params.script} \\
                -i {input.genebatch_ffn} \\
                -k {params.kmer_length} \\
                -g {params.gap_dist} \\
                -vv \\
                -o {output.prefsuffkmerfasta}
        """
