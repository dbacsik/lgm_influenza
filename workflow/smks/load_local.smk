# load_local
"""
This Snakefile takes in local data, cleans the metadata, and
filters sequences for length.
"""

# IMPORTS

# VARIABLES
out_dir = 'output/data/local/'

# RULES
rule load_local:
    input:
        fasta = 'input/data/local/Influenza_{subtype}_{segment}.fasta'
    output:
        sequences = out_dir+'{subtype}_{segment}_raw.fasta',
        metadata = out_dir+'{subtype}_{segment}_raw.tsv'
    params:
        fields = "strain segment accession originating_lab date"
    shell:
        """
        augur parse \
            --sequences {input.fasta} \
            --fields {params.fields} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fix-dates monthfirst
        """