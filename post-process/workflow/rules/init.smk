import collections
import glob
from pprint import pprint
from pathlib import Path
import os
import sys


configfile: "config.yaml"


def dir_input():
    return Path(config["input_dir"])


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
    os.path.join(f"{dir_prevoutput()}", "*/*/passing_genes/*/*_clusters.csv")
)
for g in PASSING_GENES_GLOB:
    results = g.split("/")
    PGS.append(results[-2])
    PG_GENEBATCH.append(results[-4])
    PG_BATCH.append(results[-5])


def get_passinggenes():
    return PGS


def get_passinggenes_batches():
    return PG_BATCH


def get_passinggenes_genebatches():
    return PG_GENEBATCH


def get_full_genomes():
    return DECOMP_FULL_GENOMES


def get_genomes():
    return GENOMES


def get_batches():
    return BATCHES


def fn_bakta_gff(_batch, _genome):
    return f"{dir_intermediate()}/bakta_out/{_batch}/{_genome}/{_genome}.gbff"


def fn_bakta_faa(_batch, _genome):
    return f"{dir_intermediate()}/bakta_out/{_batch}/{_genome}/{_genome}.faa"


def fn_fof(_batch):
    return f"{dir_previnput()}/{_batch}.txt"


def fn_regionfa(_batch, _genebatch, _passinggene, _contig):
    return f"{dir_intermediate()}/decompressed_regions/{_batch}/{_genebatch}/{_passinggene}/{_contig}.fasta"


def fn_genebatch_input(_genebatch):
    return f"{dir_previnput()}/{_genebatch}.ffn"


def fn_passinggenefasta(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/passinggene_fastas/{_batch}/{_genebatch}/{_passinggene}.fasta"


def fn_minimaprawout(_batch, _genebatch, _passinggene, _contig):
    return f"{dir_intermediate()}/minimap2/raw/{_batch}/{_genebatch}/{_passinggene}/{_contig}.paf"


def fn_downsampled_df(_batch, _genebatch):
    return f"{dir_intermediate()}/decompressed_genomes/{_batch}/{_genebatch}/downsampled.csv"


def fn_cluster_csv(_batch, _genebatch, _passinggene):
    return f"{dir_prevoutput()}/{_batch}/{_genebatch}/passing_genes/{_passinggene}/{_passinggene}_clusters.csv"


def fn_decompressedgenomes(_batch, _genebatch):
    return f"{dir_intermediate()}/genomes-{_batch}-{_genebatch}.txt"


def fn_decompregion_agg(_batch, _genebatch, _passinggene):
    return f"{dir_output()}/decompregions--{_batch}--{_genebatch}--{_passinggene}.txt"


def fn_decompgenome_agg(_batch, _genebatch, _passinggene):
    return f"{dir_output()}/decompgenomes--{_batch}--{_genebatch}--{_passinggene}.txt"


def fn_downsample_nwk(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/decompressed_genomes/{_batch}/{_genebatch}/{_passinggene}.nwk"


def fn_colorrange_annot(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/decompressed_genomes/{_batch}/{_genebatch}/{_passinggene}--itol/colorrange.annot"


def fn_simplebar_annot(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/decompressed_genomes/{_batch}/{_genebatch}/{_passinggene}--itol/simplebar.annot"

def fn_newtree(_batch, _genebatch, _passinggene):
    return f"{dir_intermediate()}/decompressed_genomes/{_batch}/{_genebatch}/{_passinggene}--itol/{_passinggene}--itol.nwk"

def chkpntaggregate_regionalfastas(wildcards):
    checkpoint_output = checkpoints.region_decompression.get(**wildcards).output[0]
    x = expand(
        fn_minimaprawout("{batch}", "{genebatch}", "{passinggene}", "{contig}"),
        contig=glob_wildcards(os.path.join(checkpoint_output, "{contig}.fasta")).contig,
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
