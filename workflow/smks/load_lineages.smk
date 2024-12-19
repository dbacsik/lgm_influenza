# load_local
"""
This Snakefile takes in lineage reference data.
"""

# LINEAGE DATA
"""
To have a context dataset for diverse lineages, on 2024-DEC-19, I selected at least two
sequences to represent each clade from the Nextstrain 12-year build. I selectede one
temporally early sequence and one temporally late sequence. Then, I selected a few
more sequences from Clade V1A.3a.2 to represent the diversity of that clade, since
this is all that has circulated in recent years.
"""

# IMPORTS

# VARIABLES
out_dir = 'output/build-lineages/data/'

# RULES
rule load_lineages:
    input:
        fasta = 'input/data/build-lineages/CladeReference_{subtype}_{segment}.fasta'
    output:
        sequences = out_dir+'{subtype}_{segment}_lineages.fasta',
        metadata = out_dir+'{subtype}_{segment}_lineages.tsv'
    params:
        fields = "strain segment accession originating_lab date clade"
    shell:
        """
        augur parse \
            --sequences {input.fasta} \
            --fields {params.fields} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fix-dates monthfirst
        """