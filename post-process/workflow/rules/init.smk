import collections
import glob
from pprint import pprint
from pathlib import Path
import os
import sys


configfile: "config.yaml"


def dir_intermediate():
    return Path(config["intermediate_dir"])


def dir_output():
    return config["output_dir"]


def dir_previnput():
    return Path(config["previnput_dir"])


def dir_previntermediate():
    return Path(config["previntermediate_dir"])


def dir_prevoutput():
    return Path(config["prevoutput_dir"])


PG_BATCH = []
PG_GENEBATCH = []
PGS = []
PASSING_GENES_GLOB = glob.glob(
    os.path.join(f"{dir_prevoutput()}", "*/*/passing_genes/*/*_clusters.csv.gz")
)
# print(os.path.join(f"{dir_prevoutput()}", "*/*/passing_genes/*/*_clusters.csv.gz"))
for g in PASSING_GENES_GLOB:
    # print(g)
    results = g.split("/")
    PGS.append(results[-2])
    PG_GENEBATCH.append(results[-4])
    PG_BATCH.append(results[-5])

#print(PGS)
#print(PG_GENEBATCH)
#print(PG_BATCH)


def get_passinggenes():
    return PGS

def get_newpassinggenes(wildcards):
    checkpoint_output = checkpoints.find_passinggenes.get(**wildcards).output[0]
    NEW_PASSING_GENES_GLOB = glob.glob(
        os.path.join(f"{checkpoint_output}", "passing_genes/*/*/*_clusters.csv.gz")
    )
    PG_BATCH = []
    PG_GENEBATCH = []
    PGS = []
    for g in NEW_PASSING_GENES_GLOB:
        results = g.split("/")
        PGS.append(results[-2])
        PG_GENEBATCH.append(results[-4])
        PG_BATCH.append(results[-5])
    print(PGS)
    print(PG_GENEBATCH)
    print(PG_BATCH)
    return PG_BATCH,PG_GENEBATCH,PGS




def get_passinggenes_batches():
    return PG_BATCH


def get_passinggenes_genebatches():
    return PG_GENEBATCH


def fn_fastmapraw_agg(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/insert_or_del/fastmaps--{_batch}--{_genebatch}--{_passinggene}.txt"


def fn_fastmapanalysis(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/insert_or_del/fastmap_analysis/{_batch}/{_genebatch}/{_passinggene}.fastmapanalysis.tsv"


def fn_selectedbuckets(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/insert_or_del/fastmap_analysis/{_batch}/{_genebatch}/{_passinggene}.bucket"


def fn_bakta_annot_done(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/bakta_out/{_batch}/{_genebatch}/{_passinggene}/.bakta_complete"


def fn_analyzeblocks(_batch, _genebatch, _passinggene, _bucket):
    return f"{dir_intermediate()}/insert_or_del/fastmap_analysis/{_batch}/{_genebatch}/{_passinggene}/{_bucket}.analysis.gz"


def fn_intlkmers(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/insert_or_del/internal_kmers/{_batch}/{_genebatch}/{_passinggene}_internalkmers.fa"


def fn_bwafastmap(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/insert_or_del/fastmap_raw/{_batch}/{_genebatch}/{_passinggene}"


def fn_blastidx_done(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/blast/{_batch}/{_genebatch}/{_passinggene}/.blastdb_complete"


def fn_blast_done(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/blast/{_batch}/{_genebatch}/{_passinggene}/.blast_complete"


def fn_bakta_gff(_batch, _genome):
    return f"{dir_intermediate()}/bakta_out/{_batch}/{_genome}/{_genome}.gbff"


def fn_bakta_faa(_batch, _genome):
    return f"{dir_intermediate()}/bakta_out/{_batch}/{_genome}/{_genome}.faa"


def fn_fof(_batch):
    return f"{dir_previnput()}/{_batch}.txt"


def fn_regionfadir(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/decompressed/regions/{_batch}/{_genebatch}/{_passinggene}"


def fn_regionfa(_batch, _genebatch, _passinggene, _contig):
    return fn_regionfadir(_batch, _genebatch, _passinggene) + f"{_contig}.fasta"


def fn_bucketfastmap(_batch, _genebatch, _passinggene, _bucket):
    return f"{dir_intermediate()}/insert_or_del/fastmap_analysis/{_batch}/{_genebatch}/{_passinggene}/{_bucket}.fastmap.gz"


def fn_genebatch_input(_genebatch):
    return f"{dir_previnput()}/{_genebatch}.ffn"


def fn_passinggenefasta(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/passinggene_fastas/{_batch}/{_genebatch}/{_passinggene}.fasta"


def fn_minimaprawout(_batch, _genebatch, _passinggene, _contig):
    return f"{dir_intermediate()}/minimap2/raw/{_batch}/{_genebatch}/{_passinggene}/{_contig}.paf"


def fn_blastrawout(_batch, _genebatch, _passinggene, _contig):
    return f"{dir_intermediate()}/blast/raw/{_batch}/{_genebatch}/{_passinggene}/{_contig}.blastout"


def fn_downsampled_df(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/downsampled_df/{_batch}/{_genebatch}/{_passinggene}_downsampled.csv.gz"
    


def fn_tinydownsampled_df(_batch, _genebatch, _passinggene):
    return (
        f"{dir_intermediate()}/downsampled_df/{_batch}/{_genebatch}/{_passinggene}_tinydownsampled.csv.gz",
    )


def fn_downsampled_annot_df(_batch, _genebatch, _passinggene):
    return (
        f"{dir_intermediate()}/downsampled_annot_df/{_batch}/{_genebatch}/{_passinggene}/{_passinggene}_downsampled.csv.gz",
    )


def fn_cluster_csv(_batch, _genebatch, _passinggene):
    return f"{dir_prevoutput()}/{_batch}/{_genebatch}/passing_genes/{_passinggene}/{_passinggene}_clusters.csv.gz"


def fn_clusterzero_csv(_batch, _genebatch, _passinggene):
    return (
        f"{dir_intermediate()}/downsampled_df/{_batch}/{_genebatch}/{_passinggene}/{_passinggene}_clusterzero.csv.gz",
    )


def fn_decompressedgenomes(_batch, _genebatch):
    return f"{dir_intermediate()}/genomes-{_batch}-{_genebatch}.txt"


def fn_decompregion_agg(_batch, _genebatch, _passinggene):
    return f"{dir_output()}/decompregions--{_batch}--{_genebatch}--{_passinggene}.txt"


def fn_decompgenome_agg(_batch, _genebatch, _passinggene):
    return f"{dir_output()}/decompgenomes--{_batch}--{_genebatch}--{_passinggene}.txt"


def fn_downsample_nwk(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/attotree/{_batch}/{_genebatch}/raw/{_passinggene}.nwk"


def fn_decompgenome_list(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/decompressed/genomes/large_sample/{_batch}/{_genebatch}/{_passinggene}_filelist.txt"

def fn_tinydecompgenome_list(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/decompressed/genomes/tiny_sample/{_batch}/{_genebatch}/{_passinggene}_tiny_filelist.txt"

def fn_tinydecompregion_list(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/decompressed/complete/tiny_sample/{_batch}/{_genebatch}/{_passinggene}_regiontiny_filelist.txt"

def fn_blastouts_mgedb(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/decompressed/regions/tiny_sample/{_batch}/{_genebatch}/{_passinggene}_tiny_filelist.txt"

def fn_colorrange_annot(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/attotree/{_batch}/{_genebatch}/itol/{_passinggene}--itol/colorrange.annot"


def fn_simplebar_annot(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/attotree/{_batch}/{_genebatch}/itol/{_passinggene}--itol/simplebar.annot"


def fn_itoltree(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/attotree/{_batch}/{_genebatch}/itol/{_passinggene}--itol/{_passinggene}--itol.nwk"


def fn_blastouts(_batch, _genebatch, _passinggene):
    regionfadir = fn_regionfadir(_batch, _genebatch, _passinggene)
    x = expand(
        fn_blastrawout(_batch, _genebatch, _passinggene, "{contig}"),
        contig=glob_wildcards(os.path.join(regionfadir, "{contig}.fasta")).contig,
    )
    return x


def chkpntaggregate_regionalfastas(wildcards):
    checkpoint_output = checkpoints.region_decompression.get(**wildcards).output[0]
    x = expand(
        fn_regionfa("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
        contig=glob_wildcards(os.path.join(checkpoint_output, "{contig}.fasta")).contig,
        batch=wildcards.batch,
        genebatch=wildcards.genebatch,
        passinggene=wildcards.passinggene,
    )
    return x


def chkpntaggregate_regionalfasta_dir(wildcards):
    checkpoint_output = checkpoints.genome_decompression.get(**wildcards).output[0]
    return checkpoint_output

def chkpntaggregate_passinggenes(wildcards):
    checkpoint_output = checkpoints.find_passinggenes.get(**wildcards).output[0]
    return glob.glob(os.path.join(checkpoint_output, "*/*/*_clusters.csv.gz"))


def get_filtered_gene_outputs(wildcards):
    #checkpoint_output = checkpoints.find_passinggenes.get(**wildcards).output[0]
    checkpoint_output = checkpoints.filter_passinggenes_mge.get(**wildcards).output[0]
    pattern = os.path.join(checkpoint_output, "*", "*", "*_clusters.csv.gz")
    files = glob.glob(pattern)

    outputs = []
    for f in files:
        parts = f.strip().split(os.sep)
        batch = parts[
                -3
                ]
        genebatch = parts[-2]
        basename = os.path.basename(f)
        passinggene = basename.replace("_clusters.csv.gz", "")
        outputs.append(fn_downsampled_df(batch, genebatch, passinggene))
        outputs.append(fn_downsample_nwk(batch, genebatch, passinggene))
        outputs.append(fn_itoltree(batch, genebatch, passinggene))
        outputs.append(fn_mincuts(batch, genebatch, passinggene))
        outputs.append(fn_passinggenefasta(batch, genebatch, passinggene))
        outputs.append(fn_bakta_annot_done(batch, genebatch, passinggene))
    return outputs

#def read_passing_genes(path):
    #with open(path, 'r') as o:
        

def fn_passinggenes_agg():
    return f"{dir_output()}/passinggenes.txt"

def fn_passinggene_csv(_batch, _genebatch, _passinggene):
    #return f"{dir_intermediate()}/insert_or_del/passing_genes/{_batch}/{_genebatch}/{_passinggene}_clusters.csv.gz"
    return f"{dir_intermediate()}/insert_or_del_mgedb/passing_genes/{_batch}/{_genebatch}/{_passinggene}_clusters.csv.gz"


def fn_blastouts_mgedb(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/insert_or_del_mgedb/blastouts/{_batch}/{_genebatch}/{_passinggene}_mgedb_blastouts.tsv"

def get_fastmap_output_for_bucket(wildcards):
    checkpoint_output_dir = checkpoints.run_bwafastmap.get(
        batch=wildcards.batch,
        genebatch=wildcards.genebatch,
        passinggene=wildcards.passinggene,
    ).output[
        0
    ]  # directory

    buckets = glob_wildcards(
        os.path.join(checkpoint_output_dir, "{bucket}.fastmap.gz")
    ).bucket
    return expand(
        os.path.join(checkpoint_output_dir, "{bucket}.fastmap.gz"), bucket=buckets
    )
    # path = os.path.join(checkpoint_output_dir, f"{wildcards.bucket}.fastmap.gz")
    # if not os.path.exists(path):
    # raise ValueError(
    # f"Expected fastmap file for bucket {wildcards.bucket} at {path}, but it does not exist."
    # )
    # return path


def get_fastmap_outputs(wildcards):
    checkpoint_output_dir = checkpoints.run_bwafastmap.get(
        batch=wildcards.batch,
        genebatch=wildcards.genebatch,
        passinggene=wildcards.passinggene,
    ).output[
        0
    ]  # this is the directory output

    buckets = glob_wildcards(
        os.path.join(checkpoint_output_dir, "{bucket}.fastmap.gz")
    ).bucket

    return expand(
        os.path.join(checkpoint_output_dir, "{bucket}.fastmap.gz"), bucket=buckets
    )


def analyze_bucket_outputs(wildcards):
    # Trigger re-evaluation by referencing the checkpoint
    checkpoint_dir = checkpoints.run_bwafastmap.get(
        batch=wildcards.batch,
        genebatch=wildcards.genebatch,
        passinggene=wildcards.passinggene,
    ).output[0]

    buckets = glob_wildcards(os.path.join(checkpoint_dir, "{bucket}.fastmap.gz")).bucket
    return expand(
        fn_analyzeblocks("{batch}", "{genebatch}", "{passinggene}", "{bucket}"),
        batch=wildcards.batch,
        genebatch=wildcards.genebatch,
        passinggene=wildcards.passinggene,
        bucket=buckets,
    )


def discover_buckets(wildcards):
    checkpoint_output_dir = checkpoints.run_bwafastmap.get(
        batch=wildcards.batch,
        genebatch=wildcards.genebatch,
        passinggene=wildcards.passinggene,
    ).output[0]
    return glob_wildcards(
        os.path.join(checkpoint_output_dir, "{bucket}.fastmap.gz")
    ).bucket


def chkpntaggregate_fastmap(wildcards):
    checkpoint_output = checkpoints.run_bwafastmap.get(**wildcards).output[0]
    x = expand(
        fn_bucketfastmap("{batch}", "{genebatch}", "{passinggene}", "{bucket}"),
        bucket=glob_wildcards(
            os.path.join(checkpoint_output, "{bucket}.fastmap.gz")
        ).bucket,
        batch=wildcards.batch,
        genebatch=wildcards.genebatch,
        passinggene=wildcards.passinggene,
    )
    # return "\n".join(x)
    return x


def chkpntaggregate_blastouts(wildcards):
    checkpoint_output = checkpoints.run_blast.get(**wildcards).output[0]
    x = expand(
        fn_blastrawout("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
        contig=glob_wildcards(
            os.path.join(checkpoint_output, "{contig}.blastout")
        ).contig,
        batch=wildcards.batch,
        genebatch=wildcards.genebatch,
        passinggene=wildcards.passinggene,
    )
    return x


def chkpntaggregate_genomefastas(wildcards):
    checkpoint_output = checkpoints.genome_decompression.get(**wildcards).output[0]
    x = expand(
        f"{checkpoint_output}" + "/{genome}.fa",
        genome=glob_wildcards(os.path.join(checkpoint_output, "{genome}.fa")).genome,
        batch=wildcards.batch,
        genebatch=wildcards.genebatch,
        passinggene=wildcards.passinggene,
    )
    return x


def chkpntaggregate_genomefasta_dir(wildcards):
    checkpoint_output = checkpoints.genome_decompression.get(**wildcards).output[0]
    return checkpoint_output

def chkpntaggregate_tinygenomefasta_dir(wildcards):
    checkpoint_output = checkpoints.tiny_genome_decompression.get(**wildcards).output[0]
    return checkpoint_output

def chkpntaggregate_tinyregionfasta_dir(wildcards):
    checkpoint_output = checkpoints.make_internal_fasta.get(**wildcards).output[0]
    return checkpoint_output


def fn_decompressedgenomes(_batch, _genebatch):
    return (f"{dir_intermediate()}/genomes-{_batch}-{_genebatch}.txt",)


def fn_blastidx(_batch, _genebatch, _passinggene, _contig):
    return f"{dir_intermediate()}/decompressed/regions/{_batch}/{_genebatch}/{_passinggene}/{_contig}.ndb"


def fn_selectionanalysis(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/selection/{_batch}/{_genebatch}/{_passinggene}/dnds_results.tsv"


def fn_mutationjson(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/selection/{_batch}/{_genebatch}/{_passinggene}/mutation_positions.json"


def fn_mincuts(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/minimal_cuts/{_batch}/{_genebatch}/{_passinggene}-mincuts"


def chkpntaggregate_genomes(wildcards):
    checkpoint_output = checkpoints.passinggene_cluster_decompression.get(
        **wildcards
    ).output[0]
    x = expand(
        f"{checkpoint_output}" + "/{genomes}.fa",
        genomes=glob_wildcards(os.path.join(checkpoint_output, "{genomes}.fa")).genomes,
    )
    return x


def fn_plotoutput_insertordel(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/insert_or_del/fastmap_analysis/{_batch}/{_genebatch}/{_passinggene}/all_percentage_present_by_cluster_kde.png"


def fn_percdists_insertordel(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/insert_or_del/fastmap_analysis/{_batch}/{_genebatch}/{_passinggene}/per_cluster_summary.tsv"
