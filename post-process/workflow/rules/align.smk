rule make_blastidx:
    input:
        filelist=fn_decompgenome_list("{batch}", "{genebatch}", "{passinggene}"),
        downsampled_annot_df=fn_downsampled_annot_df(
            "{batch}", "{genebatch}", "{passinggene}"
        ),
    output:
        fn_blastidx_done("{batch}", "{genebatch}", "{passinggene}"),
    conda:
        "../envs/blast.yml"
    shell:
        """
        mkdir -p $(dirname {output})
        zcat {input.downsampled_annot_df} | sed 1d | while read -r p; do
            seqname=$(echo $p | cut -d, -f2)
            inputfa=$(grep "$seqname" {input.filelist})
            makeblastdb -in $inputfa -dbtype nucl &
        done
        touch {output}
        """


rule run_blast:
    input:
        filelist=fn_decompgenome_list("{batch}", "{genebatch}", "{passinggene}"),
        downsampled_annot_df=fn_downsampled_annot_df(
            "{batch}", "{genebatch}", "{passinggene}"
        ),
        queryfasta=fn_passinggenefasta("{batch}", "{genebatch}", "{passinggene}"),
        blastidx=fn_blastidx_done("{batch}", "{genebatch}", "{passinggene}"),
    output:
        fn_blast_done("{batch}", "{genebatch}", "{passinggene}"),
    conda:
        "../envs/blast.yml"
    shell:
        """
        zcat {input.downsampled_annot_df} | sed 1d | while read -r p; do
            seqname=$(echo $p | cut -d, -f2)
            inputfa=$(grep "$seqname" {input.filelist})
            blastn -query {input.queryfasta} -subject $inputfa -outfmt 6 \\
                    > $(dirname {output})/$seqname.blastout
        done
        touch {output}
        """


# rule region_minimap:
# input:
# reffasta=fn_regionfa("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
# queryfasta=fn_passinggenefasta("{batch}", "{genebatch}", "{passinggene}"),
# output:
# fn_minimaprawout("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
# conda:
# "../envs/minimap2.yml"
# shell:
# """
# minimap2 -x splice -un -k8 {input.reffasta} {input.queryfasta} > {output}
# """
#
# rule make_blastidx:
#    input:
#        reffasta=fn_regionfa("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
#    output:
#        fn_blastidx("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
#    conda:
#        "../envs/blast.yml"
#    shell:
#        """
#        makeblastdb -in {input} -dbtype nucl
#        """
#
#
# checkpoint run_blast:
#    input:
#        reffastas=chkpntaggregate_regionalfastas,
#        #reffastadir=fn_regionfadir("{batch}", "{genebatch}", "{passinggene}"),
#        #reffasta=fn_regionfa("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
#        queryfasta=fn_passinggenefasta("{batch}", "{genebatch}", "{passinggene}"),
#    output:
#        #chkpntaggregate_blastouts,
#        #fn_blastouts("{batch}", "{genebatch}", "{passinggene}")
#        #expand(
#        #fn_blastrawout("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
#        #contig=glob_wildcards("{input.reffastadir}/{contig}.fasta").contig
#        #)
#        #expand(
#        directory(
#            f"{dir_intermediate()}" + "/blast/raw/{batch}/{genebatch}/{passinggene}"
#        ),
#        #fn_blastrawout("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
#        #contig=wildcards.contig,
#        #)
#    conda:
#        "../envs/blast.yml"
#    shell:
#        #"""
#        #blastn -query {input.queryfasta} -subject {input.reffasta} -outfmt 6 > {output}
#        #"""
#        """
#        x="{input.reffastas}"
#        mkdir -p {output}
#        for i in $x; do
#            blastn -query {input.queryfasta} -subject $i -outfmt 6 \\
#                    > {output}/$(basename $i .fasta).blastout
#        done
#        """
