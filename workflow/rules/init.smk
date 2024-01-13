# from snakemake.utils import validate
import collections
import glob
from pprint import pprint
from pathlib import Path
import os


configfile: "config.yaml"


##### load config and sample sheets #####


def dir_input():
    return Path(config["input_dir"])


def dir_intermediate():
    return Path(config["intermediate_dir"])


def dir_output():
    return config["output_dir"]


# extract sample name from a path


def _get_sample_from_fn(x):
    suffixes = ["fa", "fasta", "fna", "ffa"]

    b = os.path.basename(x)
    sample, _, suffix = b.rpartition(".")
    assert suffix in suffixes, f"Unknown suffix of source files ({suffix} in {x})"
    return sample


# compute main dict for batches
# TODO: if executed in cluster mode, every job will recompute this BATCHES_FN variable when submitted
# TODO: this is because in cluster mode each job is ran as <get the actual snakemake command line if needed>
# TODO: this makes us include this file and thus recompute BATCHES_FN
# TODO: might be a good idea to serialise BATCHES_FN to disk and read from it, instead of recomputing it every time
# TODO: it might hammer the disk in cluster envs, depending on the number of batches
BATCHES_FN = {}
SAMPLE_TO_BATCHES = {}
res = dir_input().glob("*.txt")
for x in res:
    b = os.path.basename(x)
    if not b.endswith(".txt"):
        continue
    batch = b[:-4]

    BATCHES_FN[batch] = {}
    with open(x) as f:
        for y in f:
            sample_fn = y.strip()
            if not sample_fn.endswith(".tar.xz"):
                continue
            sample_base = os.path.basename(sample_fn).split(".tar.xz")[0]
            SAMPLE_TO_BATCHES[sample_fn] = batch
            BATCHES_FN[batch][sample_base] = sample_fn

GENEBATCHES_FN = []
res = dir_input().glob("*.ffn")
for x in res:
    b = os.path.basename(x)
    if not b.endswith(".ffn"):
        continue
    batch = b[:-4]

    GENEBATCHES_FN.append(batch)


assert (
    len(BATCHES_FN) != 0
), f"\nERROR: No input genomes provided. Please provide at least one batch in '{dir_input()}/'. Must have the '.txt' extension.\n"

assert (
    len(GENEBATCHES_FN) != 0
), f"\nERROR: No input genes provided. Please provide at least one batch in '{dir_input()}/'. Must have the '.ffn' extension.\n"


## BATCHES


def get_genome_batches():
    return BATCHES_FN.keys()


def get_samples():
    return SAMPLE_TO_BATCHES.keys()


def get_batch(_sample):
    return SAMPLE_TO_BATCHES[_sample]


def get_sample_paths(_batch, _sample):
    return BATCHES_FN[_batch][_sample]


def get_gene_batches():
    return GENEBATCHES_FN


#####################################
# Global files for individual batches
#####################################


def fn_batchidxdir(_batch, _sample):
    return f"{dir_intermediate()}/bwaidx/{_batch}/{_sample}"


def aggregate_fastmap_dists(wildcards):
    checkpoint_output = checkpoints.make_bwaidx.get(**wildcards).output[0]
    return expand(
        f"{dir_intermediate()}"
        + "/fastmap/distances/{batch}/{genebatch}/{genebatch}-{bucket}-distances.tsv",
        batch=wildcards.batch,
        genebatch=get_gene_batches(),
        bucket=glob_wildcards(os.path.join(checkpoint_output, "{bucket}.bwt")).bucket,
    )


def aggregate_filter_dists(wildcards):
    checkpoint_output = checkpoints.make_bwaidx.get(**wildcards).output[0]
    return expand(
        f"{dir_intermediate()}"
        + "/kmerdists/{batch}/{genebatch}/{bucket}-filterdists.csv",
        batch=wildcards.batch,
        genebatch=get_gene_batches(),
        bucket=glob_wildcards(os.path.join(checkpoint_output, "{bucket}.bwt")).bucket,
    )


def aggregate_filter_distdirs(wildcards):
    checkpoint_output = checkpoints.make_bwaidx.get(**wildcards).output[0]
    return expand(
        f"{dir_intermediate()}" + "/kmerdists/{batch}/{genebatch}",
        batch=wildcards.batch,
        genebatch=get_gene_batches(),
    )


def aggregate_passing_genes(wildcards):
    # def aggregate_passing_genes():
    checkpoint_output = checkpoints.cluster_dists.get(**wildcards).output[0]
    # checkpoint_output = checkpoints.cluster_dists.get().output[0]
    return expand(
        f"{checkpoint_output}" + "/{passing_genes}_clusters.csv",
        passing_genes=glob_wildcards(
            os.path.join(checkpoint_output, "{passing_genes}_clusters.csv")
        ).passing_genes,
    )


def fn_prefsuffkmer(_genebatch):
    return f"{dir_intermediate()}/prefsuffkmers/{_genebatch}.ffn"


def fn_bwaidxbucket(_batch, _bucket):
    return f"{dir_intermediate()}/bwa_indices/{_batch}/{_bucket}.bwt"


def fn_fastmapraw(_batch, _bucket, _genebatch):
    return f"{dir_intermediate()}/fastmap/raw/{_batch}/{_genebatch}/{_genebatch}-{_bucket}-raw.fastmap"


def fn_fastmapprocess(_batch, _bucket, _genebatch):
    return f"{dir_intermediate()}/fastmap/processed/{_batch}/{_genebatch}/{_genebatch}-{_bucket}-processed.fastmap"


def fn_fastmapdists(_batch, _bucket, _genebatch):
    return f"{dir_intermediate()}/fastmap/distances/{_batch}/{_genebatch}/{_genebatch}-{_bucket}-distances.tsv"


def fn_parsedists(_batch, _bucket, _genebatch):
    return f"{dir_intermediate()}/kmerdists/{_batch}/{_genebatch}/{_bucket}-filterdists.csv"


def fn_downsampled_df(_batch, _genebatch, _passinggene):
    passing_gene = os.path.basename(_passinggene).split("_clusters.csv")[0]
    return f"{dir_intermediate()}/decompressed_genomes/{_batch}/{_genebatch}/{passing_gene}-downsampled.csv"


def fn_listoutputteddfs(passing_file_list):
    y = []
    for passing_file in passing_file_list:
        genes = []
        batch = []
        genebatch = []
        with open(passing_file, "r") as f:
            init_list = f.read().split()
            for i in init_list:
                parts = i.split("/")
                genes.append(parts[-2])
                genebatch.append(parts[-4])
                batch.append(parts[-5])
        for i in range(len(genes)):
            y.append(
                f"{dir_intermediate()}/decompressed_genomes/{batch[i]}/{genebatch[i]}/{genes[i]}-downsampled.csv"
            )
    return y


def fn_prefsuffkmer(_genebatch):
    return f"{dir_intermediate()}/prefsuffkmers/{_genebatch}-prefsuff-k{config['kmer_length']}-g{config['gap_distance']}.ffn"


## WILDCARD FUNCTIONS


# get source file path
def w_sample_source(wildcards):
    batch = wildcards["batch"]
    sample = wildcards["sample"]
    fn = BATCHES_FN[batch][sample]
    return fn


## OTHER FUNCTIONS


# generate file list from a list of identifiers (e.g., leaf names -> assemblies names)
def generate_file_list(input_list_fn, output_list_fn, filename_function):
    with open(input_list_fn) as f:
        with open(output_list_fn, "w") as g:
            for x in f:
                x = x.strip()
                fn0 = filename_function(x)  # top-level path
                fn = os.path.relpath(fn0, os.path.dirname(output_list_fn))
                g.write(fn + "\n")


def load_list(fn):
    try:
        with open(fn) as f:
            return [x.strip() for x in f]
    except FileNotFoundError:
        print(f"File not found {fn}, using empty list")
        return []
