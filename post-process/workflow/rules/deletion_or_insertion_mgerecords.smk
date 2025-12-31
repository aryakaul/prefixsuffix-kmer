checkpoint make_internal_fasta:
    output:
        outdir=directory(
            f"{dir_intermediate()}"
            + "/decompressed/complete/downsample/{batch}/{genebatch}/{passinggene}"
        ),
    input:
        downsampled_df = fn_downsampled_df("{batch}", "{genebatch}", "{passinggene}"),
        fof=fn_fof("{batch}"),
    params:
        script = Path(workflow.basedir) / "scripts/decompression",
        internal_pad = config["internal_pad"],
        KMERSIZE = config['kmer_length'],
        FLANKKB = ['flankkb_decompress']
    conda:
        "../envs/biopython.yml"
    shell:
        """
        {params.script} \
            -i {input.downsampled_df} \
            -f {input.fof} \
            -e {output.outdir} \
            -x {params.FLANKKB} \
            -k {params.KMERSIZE} \
            -g {params.internal_pad} \
            -vvv
        """

rule build_tinyfilelist_region:
    output:
        filelist=fn_tinydecompregion_list("{batch}", "{genebatch}", "{passinggene}"),
    input:
        chkpntaggregate_tinyregionfasta_dir,
    shell:
        r"""
        set -euo pipefail
        shopt -s nullglob

        fasta_files=( {input}/*/*.removed_regions.fasta )
        # If none, that's an error
        if [ ${{#fasta_files[@]}} -eq 0 ]; then
            echo "ERROR: no .fasta files found in {input}" >&2
            exit 1
        fi

        # Otherwise write the list
        printf "%s\n" "${{fasta_files[@]}}" > {output.filelist}
        """

rule blast_mgedb:
    output:
        blastout=fn_blastouts_mgedb("{batch}", "{genebatch}", "{passinggene}"),
    input:
        filelist=fn_tinydecompregion_list("{batch}", "{genebatch}", "{passinggene}"),
    params:
        mgedb=Path(workflow.basedir)
        / "../supp_data/mgedb_prophagedb/mgedb_prophagedb.fna",
    conda:
        "../envs/blast.yml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.blastout})
        : > {output.blastout}

        # Add qlen to the header
        echo -e "genome\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen" \
            > {output.blastout}

        while read fasta; do
            genome=$(basename "$fasta" .fasta)
            # Single-record FASTA assumed; compute its length for the no-hit case
            qlen=$(awk 'BEGIN{{L=0}} /^>/ {{next}} {{L+=length($0)}} END{{print L}}' "$fasta")

            # Include qlen in outfmt (one value per hit)
            hits=$(blastn -query "$fasta" -db {params.mgedb} \
                     -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen")

            if [ -n "$hits" ]; then
                # Note doubled braces to escape awk braces for Snakemake
                echo "$hits" | awk -v g="$genome" '{{print g"\t"$0}}' >> {output.blastout}
            else
                # Emit a no-hit row with qlen so denominator is correct
                echo -e "${{genome}}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t${{qlen}}" >> {output.blastout}
            fi
        done < {input.filelist}
        """

checkpoint filter_passinggenes_mge:
    output:
        directory(f"{dir_intermediate()}/insert_or_del_mgedb/passing_genes"),
    input:
        blastout = expand(
            fn_blastouts_mgedb("{batch}", "{genebatch}", "{passinggene}"),
            zip,
            batch=get_passinggenes_batches(),
            genebatch=get_passinggenes_genebatches(),
            passinggene=get_passinggenes(),
        ),
        downsampled_df = expand(
            fn_tinydownsampled_df("{batch}", "{genebatch}", "{passinggene}"),
            zip,
            batch=get_passinggenes_batches(),
            genebatch=get_passinggenes_genebatches(),
            passinggene=get_passinggenes(),
        ),
    params:
        intermediate=f"{dir_intermediate()}",
        prevoutput=f"{dir_prevoutput()}",
        max_frac=config["threshold"],
        script=Path(workflow.basedir) / "scripts/filter_mge_dominated",
    conda:
        "../envs/pandas.yml"
    shell:
        """
        mkdir -p {output}
        for i in {input.blastout}; do

            batch=$(basename $(dirname $(dirname $i)))
            genebatch=$(basename $(dirname $i))
            gene=$(basename $i _mgedb_blastouts.tsv)

            tinydf={params.intermediate}/downsampled_df/"$batch"/"$genebatch"/"$gene"_tinydownsampled.csv.gz


            outdir={output}/$batch/$genebatch
            mkdir -p $outdir

            filtered=$outdir/"$gene"_filtered.tsv

            # Run Python filter
            {params.script} --blast $i --csvfile $tinydf --output "$filtered" --max_frac {params.max_frac}

            fail=$(
              awk -F'\t' 'NR==1 {{next}} {{
                  c=$2;                        # Cluster column
                  tot[c]++
                  if (tolower($5)=="true") dom[c]++   # dominated column
                }}
                END {{
                  for (c in tot) {{
                    if (dom[c]*2 > tot[c]) {{   # strict majority
                      print 1; exit            # FAIL
                    }}
                  }}
                  print 0                      # PASS
                }}' "$filtered"
            )
            if [ "$fail" -eq 1 ]; then
                echo "FAIL: $gene excluded due to MGE dominance" >&2
            else
                logger "PASS: $gene not dominated by MGE hits"
                ln -s {params.prevoutput}/$batch/$genebatch/passing_genes/$gene/"$gene"_clusters.csv.gz \
                    $outdir
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

