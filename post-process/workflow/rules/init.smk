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

# print(PGS)


def get_passinggenes():
    return PGS


def get_passinggenes_batches():
    return PG_BATCH


def get_passinggenes_genebatches():
    return PG_GENEBATCH


def fn_bakta_annot_done(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/bakta_out/{_batch}/{_genebatch}/{_passinggene}/.bakta_complete"


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
    return f"{dir_intermediate()}/decompressed_regions/{_batch}/{_genebatch}/{_passinggene}"


def fn_regionfa(_batch, _genebatch, _passinggene, _contig):
    return f"{dir_intermediate()}/decompressed_regions/{_batch}/{_genebatch}/{_passinggene}/{_contig}.fasta"


def fn_genebatch_input(_genebatch):
    return f"{dir_previnput()}/{_genebatch}.ffn"


def fn_passinggenefasta(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/passinggene_fastas/{_batch}/{_genebatch}/{_passinggene}.fasta"


def fn_minimaprawout(_batch, _genebatch, _passinggene, _contig):
    return f"{dir_intermediate()}/minimap2/raw/{_batch}/{_genebatch}/{_passinggene}/{_contig}.paf"


def fn_blastrawout(_batch, _genebatch, _passinggene, _contig):
    return f"{dir_intermediate()}/blast/raw/{_batch}/{_genebatch}/{_passinggene}/{_contig}.blastout"


def fn_downsampled_df(_batch, _genebatch, _passinggene):
    return (
        f"{dir_intermediate()}/downsampled_df/{_batch}/{_genebatch}/{_passinggene}_downsampled.csv.gz",
    )


def fn_downsampled_annot_df(_batch, _genebatch, _passinggene):
    return (
        f"{dir_intermediate()}/downsampled_annot_df/{_batch}/{_genebatch}/{_passinggene}/{_passinggene}_downsampled.csv.gz",
    )


def fn_cluster_csv(_batch, _genebatch, _passinggene):
    return f"{dir_prevoutput()}/{_batch}/{_genebatch}/passing_genes/{_passinggene}/{_passinggene}_clusters.csv.gz"


def fn_decompressedgenomes(_batch, _genebatch):
    return f"{dir_intermediate()}/genomes-{_batch}-{_genebatch}.txt"


def fn_decompregion_agg(_batch, _genebatch, _passinggene):
    return f"{dir_output()}/decompregions--{_batch}--{_genebatch}--{_passinggene}.txt"


def fn_decompgenome_agg(_batch, _genebatch, _passinggene):
    return f"{dir_output()}/decompgenomes--{_batch}--{_genebatch}--{_passinggene}.txt"


def fn_downsample_nwk(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/attotree/{_batch}/{_genebatch}/raw/{_passinggene}.nwk"


def fn_decompgenome_list(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/decompressed_genomes/{_batch}/{_genebatch}/{_passinggene}_filelist.txt"


def fn_colorrange_annot(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/attotree/{_batch}/{_genebatch}/itol/{_passinggene}--itol/colorrange.annot"


def fn_simplebar_annot(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/attotree/{_batch}/{_genebatch}/itol/{_passinggene}--itol/simplebar.annot"


def fn_itoltree(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/attotree/{_batch}/{_genebatch}/itol/{_passinggene}--itol/{_passinggene}--itol.nwk"


def fn_blastouts(_batch, _genebatch, _passinggene):
    regionfadir = fn_regionfadir(_batch, _genebatch, _passinggene)
    print(regionfadir)
    x = expand(
        fn_blastrawout(_batch, _genebatch, _passinggene, "{contig}"),
        contig=glob_wildcards(os.path.join(regionfadir, "{contig}.fasta")).contig,
    )
    return x


def chkpntaggregate_regionalfastas(wildcards):
    checkpoint_output = checkpoints.region_decompression.get(**wildcards).output[0]
    x = expand(
        # fn_minimaprawout("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
        fn_regionfa("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
        contig=glob_wildcards(os.path.join(checkpoint_output, "{contig}.fasta")).contig,
        batch=wildcards.batch,
        genebatch=wildcards.genebatch,
        passinggene=wildcards.passinggene,
    )
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


def fn_decompressedgenomes(_batch, _genebatch):
    return (f"{dir_intermediate()}/genomes-{_batch}-{_genebatch}.txt",)


def fn_blastidx(_batch, _genebatch, _passinggene, _contig):
    return f"{dir_intermediate()}/decompressed_regions/{_batch}/{_genebatch}/{_passinggene}/{_contig}.ndb"


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
