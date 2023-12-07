rule makebwaidx:
    output:
        bwaindex=fn_idxpersample(_batch="{batch}", _sample="{sample}"),
    input:
        fasta=w_sample_source,
    conda:
        "../envs/bwamem.yml"
    shell:
        """
        mkdir -p $(dirname {output})
        bwa index -p $(dirname {output})/$(basename {output} .pac) {input.fasta}
        """


rule fastmap:
    output:
        bwafastmap=fn_fastmapraw(_sample="{sample}", _batch="{batch}", _genebatch="{genebatch}"),
    input:
        prefsuffkmers = fn_prefsuffkmer(_genebatch="{genebatch}"),
        bwaindex=fn_idxpersample(_batch="{batch}", _sample="{sample}"),
    params:
        kmer_length=config["kmer_length"],
    conda:
        "../envs/bwamem.yml"
    shell:
        """
        idx=$(dirname {input.bwaindex})/$(basename {input.bwaindex} .pac)
        bwa fastmap -w 99999 -l {params.kmer_length} $idx {input.prefsuffkmers} > {output}
        """
