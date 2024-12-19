# LGM_INFLUENZA
"""
This Snakefile organizes the LGM Influenza project. It uses Nextstrain for many steps 
and it is based on the one in the 'quickstart' repository from the public avian-flu 
build: https://github.com/nextstrain/avian-flu/blob/master/quickstart-build/Snakefile
"""

# IMPORTS

# VARIABLES
SUBTYPES = ['B-Victoria']
SEGMENTS = ['HA', 'NA']

# ALL RULE
rule all:
    input:
        ingest_local = expand(
            'output/data/local/{subtype}_{segment}_filterLength.tsv',
            subtype=SUBTYPES, segment=SEGMENTS),

# Includes
include: "workflow/smks/load_local.smk"