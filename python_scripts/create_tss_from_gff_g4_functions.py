import pandas as pd
import re
import matplotlib 
import gzip
import os
import pybedtools
import gff_g4_functions

base_dir="/uufs/chpc.utah.edu/common/home/u0674686/ncbi_genomes/ncbi-genomes-representive1200/"
gff_g4_functions.create_tss_from_gff3(base_dir)
