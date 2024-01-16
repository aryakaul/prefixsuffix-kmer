import collections
import glob
from pprint import pprint
from pathlib import Path
import os
import sys


configfile: "config.yaml"

#def dir_input():
    #return Path(config["input_dir"])

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

DECOMP_FULL_GENOMES = dir_prevoutput().glob("decompressed_genomes/genomes/*/*.fa")
BATCHES = []
GENOMES = []
for i in DECOMP_FULL_GENOMES:
    BATCHES.append(os.path.basename(os.path.dirname(i)))
    GENOMES.append(os.path.basename(i).split('.fa')[0])

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

def fn_downsampled_df(_batch, _genebatch):
    return f"{dir_intermediate()}/decompressed_genomes/{_batch}/{_genebatch}/downsampled.csv"

def fn_decompressedgenomes(_batch, _genebatch):
    return f"{dir_intermediate()}/genomes-{_batch}-{_genebatch}.txt",

def chkpntaggregate_genomes(wildcards):
    checkpoint_output = checkpoints.passinggene_cluster_decompression.get(**wildcards).output[0]
    x = expand(
        f"{checkpoint_output}" + "/{genomes}.fa",
        genomes = glob_wildcards(
            os.path.join(checkpoint_output, "{genomes}.fa")
        ).genomes
    )
    return x


