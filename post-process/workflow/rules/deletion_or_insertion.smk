rule extract_internal_kmers:
    output:
        internal_kmers=fn_intlkmers("{batch}", "{genebatch}", "{passinggene}"),
    input:
        genomefastadir=chkpntaggregate_tinygenomefasta_dir,
        filelist=fn_tinydecompgenome_list("{batch}", "{genebatch}", "{passinggene}"),
        downsampled_df=fn_tinydownsampled_df("{batch}", "{genebatch}", "{passinggene}"),
    params:
        num_internal=config["number_annotate"],
        script=Path(workflow.basedir) / "scripts/get_internal_kmers",
        KMERSIZE=config["kmer_length"],
        fastadir=Path(
            f"{dir_intermediate()}"
            + "/decompressed_genomes/{batch}/{genebatch}/{passinggene}"
        ),
        samples_per_cluster=config["num_filter_genomes"],
        kmers_per_genome=30,
    threads: 8
    conda:
        "../envs/biopython.yml"
    shell:
        """
        mkdir -p $(dirname {output})
        {params.script} \\
            --csv {input.downsampled_df} \\
            --fasta_dir {input.genomefastadir} \\
            --samples_per_cluster {params.samples_per_cluster} \\
            --kmers_per_genome {params.kmers_per_genome} \\
            --kmer_length {params.KMERSIZE} \\
            -j {threads} \\
            -v \\
            --output {output.internal_kmers}
        """


rule select_buckets:
    input:
        clustercsv=fn_cluster_csv("{batch}", "{genebatch}", "{passinggene}"),
    output:
        buckets=fn_selectedbuckets("{batch}", "{genebatch}", "{passinggene}"),
    params:
        top_n=30,
        script=Path(workflow.basedir) / "scripts/select_buckets",
    conda:
        "../envs/biopython.yml"
    shell:
        """
        mkdir -p $(dirname {output})
        {params.script} \\
            {input} \\
            {params.top_n} \\
            {output}
        """


checkpoint run_bwafastmap:
    output:
        outdir=directory(
            f"{dir_intermediate()}"
            + "/insert_or_del/fastmap_raw/{batch}/{genebatch}/{passinggene}"
        ),
    input:
        internal_kmers=fn_intlkmers("{batch}", "{genebatch}", "{passinggene}"),
        buckets=fn_selectedbuckets("{batch}", "{genebatch}", "{passinggene}"),
    params:
        KMERSIZE=config["kmer_length"],
        indices=config["fastmap_indices"],
    conda:
        "../envs/bwamem.yml"
    threads: 4
    shell:
        """
        mkdir -p {output.outdir}
        
        #cat {input.buckets} | parallel --jobs {threads} --halt soon,fail=1 \
            #'bwa fastmap -w 99999 -l {params.KMERSIZE} {params.indices}/{{}} {input.internal_kmers} | gzip > {output.outdir}/{{}}.fastmap.gz'
        while IFS= read -r bucket; do
            bwa fastmap -w 99999 -l {params.KMERSIZE} {params.indices}/"$bucket" {input.internal_kmers} | gzip > {output.outdir}/"$bucket".fastmap.gz
        done < {input.buckets}
        #wait
        """

checkpoint run_bwafastmap_isdb:
    output:
        outdir=directory(
            f"{dir_intermediate()}"
            + "/insert_or_del_isdb/fastmap_raw/{batch}/{genebatch}/{passinggene}"
        ),
    input:
        internal_kmers=fn_intlkmers("{batch}", "{genebatch}", "{passinggene}"),
    params:
        KMERSIZE=config["kmer_length"],
        indices=config["fastmap_indices"],
    conda:
        "../envs/bwamem.yml"
    threads: 4
    shell:
        """
        mkdir -p {output.outdir}
        while IFS= read -r bucket; do
            bwa fastmap -w 99999 -l {params.KMERSIZE} {params.indices}/"$bucket" {input.internal_kmers} | gzip > {output.outdir}/"$bucket".fastmap.gz
        done < {input.buckets}
        """

rule aggregate_fastmapout:
    output:
        fn_fastmapraw_agg("{batch}", "{genebatch}", "{passinggene}"),
    input:
        chkpntaggregate_fastmap,
    shell:
        """
        echo {input} > {output}
        """


rule analyze_bwafastmap:
    output:
        fastmap_analysis=fn_fastmapanalysis("{batch}", "{genebatch}", "{passinggene}"),
    input:
        internal_kmers=fn_intlkmers("{batch}", "{genebatch}", "{passinggene}"),
        clustercsv=fn_cluster_csv("{batch}", "{genebatch}", "{passinggene}"),
        fastmap=rules.run_bwafastmap.output,
    params:
        KMERSIZE=config["kmer_length"],
        script=Path(workflow.basedir) / "scripts/internal_kmers_bwaout",
        presence_threshold=0.3,
    threads: 16
    conda:
        "../envs/biopython.yml"
    shell:
        """
        {params.script} \\
            --query_fasta {input.internal_kmers} \\
            --fastmap_dir {input.fastmap} \\
            --output_tsv {output.fastmap_analysis} \\
            --presence_threshold {params.presence_threshold} \\
            --num_processes {threads} \\
            -v
        """


rule find_percentage_dist:
    output:
        #plotfn=fn_plotoutput_insertordel("{batch}", "{genebatch}", "{passinggene}"),
        percentagedists=fn_percdists_insertordel("{batch}", "{genebatch}", "{passinggene}"),
    input:
        fastmap_tsv=fn_fastmapanalysis("{batch}", "{genebatch}", "{passinggene}"),
        passingcsv=fn_cluster_csv("{batch}", "{genebatch}", "{passinggene}"),
        bucket=fn_selectedbuckets("{batch}", "{genebatch}", "{passinggene}"),
    params:
        threshold=config["threshold"],
        script=Path(workflow.basedir) / "scripts/find_cluster_summary_insertordel",
    conda:
        "../envs/pandas.yml"
    shell:
        """

        OUTPUTDIR=$(dirname {output.percentagedists})
        {params.script} \\
            {input.fastmap_tsv} \\
            {input.passingcsv} \\
            {input.bucket} \\
            --plot_output $OUTPUTDIR \\
            -v
        """

checkpoint find_passinggenes:
    output:
        directory(
            f"{dir_intermediate()}"
            + "/insert_or_del/passing_genes"
        ),
    input:
        expand(
            fn_percdists_insertordel("{batch}", "{genebatch}", "{passinggene}"),
            zip,
            batch=get_passinggenes_batches(),
            genebatch=get_passinggenes_genebatches(),
            passinggene=get_passinggenes(),
        ),
    params:
        threshold=config["threshold"],
        prevoutput=f"{dir_prevoutput()}"
    conda:
        "../envs/pandas.yml"
    shell:
        """
        mkdir -p {output}
        for i in {input}; do
            genebatch=$(basename $(dirname $(dirname $i)))
            batch=$(basename $(dirname $(dirname $(dirname $i))))
            mkdir -p {output}/$batch/$genebatch
            gene=$(basename $(dirname $i))
            output_csv={params.prevoutput}/$batch/$genebatch/passing_genes/$gene/"$gene"_clusters.csv.gz

            if [[ ! -f "$i" ]]; then
                echo "Missing summary file: $i" >&2
                continue
            fi

            fail=0
            threshold={params.threshold}
            while IFS=$'\t' read -r cluster numgenomes mean frac; do
                if [[ "$cluster" != "Cluster" && "$cluster" != "0" && "$cluster" != "-1" ]]; then
                    above=$(echo "$frac > $threshold" | bc -l)
                    if (( above )); then
                        fail=1
                        break
                    fi
                fi
            done < "$i"

            if (( fail == 0 )); then
                ln -s "$output_csv" {output}/$batch/$genebatch 
            else
                echo "FAIL: $gene excluded due to high overlap in cluster separation" >&2
            fi
        done
        """

checkpoint find_passinggenes_isdb:
    output:
        directory(
            f"{dir_intermediate()}"
            + "/insert_or_del_isdb/passing_genes"
        ),
    input:
        expand(
            fn_percdists_insertordel("{batch}", "{genebatch}", "{passinggene}"),
            zip,
            batch=get_passinggenes_batches(),
            genebatch=get_passinggenes_genebatches(),
            passinggene=get_passinggenes(),
        ),
    params:
        threshold=config["threshold"],
        prevoutput=f"{dir_prevoutput()}"
    conda:
        "../envs/pandas.yml"
    shell:
        """
        mkdir -p {output}
        for i in {input}; do
            genebatch=$(basename $(dirname $(dirname $i)))
            batch=$(basename $(dirname $(dirname $(dirname $i))))
            mkdir -p {output}/$batch/$genebatch
            gene=$(basename $(dirname $i))
            output_csv={params.prevoutput}/$batch/$genebatch/passing_genes/$gene/"$gene"_clusters.csv.gz

            if [[ ! -f "$i" ]]; then
                echo "Missing summary file: $i" >&2
                continue
            fi

            fail=0
            threshold={params.threshold}
            while IFS=$'\t' read -r cluster numgenomes mean frac; do
                if [[ "$cluster" != "Cluster" && "$cluster" != "0" && "$cluster" != "-1" ]]; then
                    above=$(echo "$frac > $threshold" | bc -l)
                    if (( above )); then
                        fail=1
                        break
                    fi
                fi
            done < "$i"

            if (( fail == 0 )); then
                ln -s "$output_csv" {output}/$batch/$genebatch 
            else
                echo "FAIL: $gene excluded due to high overlap in cluster separation" >&2
            fi
        done
        """



rule aggregate_downsample_df:
    input:
        get_filtered_gene_outputs

rule collect_passinggenes:
    output:
        fn_passinggenes_agg()
    input:
        chkpntaggregate_passinggenes
    shell:
        """
        echo {input} > {output}
        """
