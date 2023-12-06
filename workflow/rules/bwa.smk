rule makebwaidx:
    output:
        bwaindex=fn_idxpersample(_sample="{sample}"),
    input:
        sample=fn_samplefasta(_batch="{batch}", _sample="{sample}"),
    conda:
        "../envs/bwamem.yml"
    params:
        kmer_length=config['kmer_length']
    shell:
        """
        mkdir -p $(dirname {output})
        bwa index -p $(basename {output} .pac) -k {params.kmer_length} {sample}
        """

rule fastmap:
    output:
        bwafastmap=fn_fastmapraw(_sample="{sample}"),
    input:
        prefsuffkmers=fn_prefsuffkmer(_genebatch="{genebatch}"),
        sample=fn_samplefasta(_batch="{batch}", _sample="{sample}"),
    params:
        kmer_length=config["kmer_length"],
    conda:
        "../envs/bwamem.yml"
    shell:
        """
        bwa fastmap -w 99999 -l {params.kmer_length} {input.sample} {input.prefsuffkmers} > {output}
        """
