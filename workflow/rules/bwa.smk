rule makebwaidx:
    output:
        bwaindex=fn_idxpersample(_sample="{sample}"),
    input:
        sample=fn_sample(_sample="{sample}"),
    conda:
        "../envs/bwamem.yml"
    shell:
        """
        bwa index {input}
        """

rule fastmap:
    output:
        bwafastmap=fn_fastmapraw(_sample="{sample}"),
    input:
        prefsuffkmers=fn_prefsuffkmer(_gene="{gene}"),
        sample=fn_sample(_sample="{sample}"),
    params:
        kmer_length=config["kmer_length"],
    conda:
        "../envs/bwamem.yml"
    shell:
        """
        bwa fastmap -w 99999 -l {params.kmer_length} {input.sample} {input.prefsuffkmers} > {output}
        """
