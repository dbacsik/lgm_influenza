# load_local
"""
This Snakefile takes in local data, cleans the metadata, and
filters sequences for length.
"""

# IMPORTS

# VARIABLES
workdir: '/workspaces/flu_pisaac-2024/'
out_dir = 'output/data/local/'
min_length = { # Minimum length for each segment, copied from Nextstrain
    'HA': 1500,
    'NA': 1400
}

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

rule filter_length:
    input:
        sequences = rules.load_local.output.sequences,
        metadata = rules.load_local.output.metadata
    output:
        sequences = out_dir+'{subtype}_{segment}_filterLength.fasta',
        metadata = out_dir+'{subtype}_{segment}_filterLength.tsv'
    params:
        min_length = lambda wc: min_length[wc.segment]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --min-length {params.min_length} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}
        """

rule parse_country:
    input:
        metadata = rules.filter_length.output.metadata
    output:
        metadata = out_dir+'{subtype}_{segment}_parseCountry.tsv'
    notebook:
        '../notebooks/parse_country.py.ipynb'

rule add_local_flag:
    input:
        metadata = rules.parse_country.output.metadata
    output:
        metadata = out_dir+'{subtype}_{segment}_localFlag.tsv'
    shell:
        """
        awk 'BEGIN {{FS=OFS="\t"}} {{print $0, "local"}}' \
            {input.metadata} > {output.metadata}
        """