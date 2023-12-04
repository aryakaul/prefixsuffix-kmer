from snakemake.utils import validate
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
    if b.endswith(".gz"):
        b = b[:-3]
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
            if sample_fn:
                sample = _get_sample_from_fn(sample_fn)
                BATCHES_FN[batch][sample] = sample_fn

assert (
    len(BATCHES_FN) != 0
), f"\nERROR: No input files provided. Please provide at least one batch in '{dir_input()}/'.\n"


## BATCHES


def get_batches():
    return BATCHES_FN.keys()


#####################################
# Global files for individual batches
#####################################


def fn_tree_clean(_batch):
    return f"{dir_intermediate()}/tree/{_batch}.nwk"


def fn_tree_sorted(_batch):
    return f"{dir_intermediate()}/tree/{_batch}.nwk"


def fn_tree_dirty(_batch):
    return f"{dir_intermediate()}/tree/{_batch}.nwk_dirty"


def fn_leaves_sorted(_batch):
    return f"{dir_intermediate()}/tree/{_batch}.leaves"

def fn_nodes_sorted(_batch):
    return f"{dir_intermediate()}/tree/{_batch}.nodes"

#####################################
# Files for individual samples      #
#####################################

def fn_idxpersample(_sample):
    pass

def fn_sample(_sample):
    pass

def fn_fastmapraw(_sample):
    pass

def fn_prefsuffkmer(_gene):
    pass

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
