# build_lineages
"""
This Snakefile assigns lineages to our influenza samples. It builds a minimal tree
of representative sequences annotated with each lineage and then infers
the lineage for each new sample.
"""

# IMPORTS
include: 'load_local.smk'
include: 'load_lineages.smk'

# VARIABLES
build_dir = 'output/build-lineages/'

ref_files = {
    'B-Victoria': {
        'HA': {
            'FASTA': 'input/reference/OQ202354_B-Victoria_HA.fasta',
            'gbk': 'input/reference/OQ202354_B-Victoria_HA.gb',
            'clade_muts': 'input/data/build-lineages/Influenza_B-Victoria_HA_clades.tsv'
        },
        'NA': {
            'FASTA': 'input/reference/OQ202364_B-Victoria_NA.fasta',
            'gbk': 'input/reference/OQ202364_B-Victoria_NA.gb'
        }
    }
}

# RULES
# DATA PREPARATION
## COMBINE DATA

"""This rule concatenates the local sequences with the reference dataset."""
rule cat_fastas:
    message:
        """
        Concatenating the lineage reference FASTA and the local FASTA.\n
        Reference FASTA: {input.reference}\n
        Local FASTA: {input.local}
        """
    input:
        reference = rules.load_lineages.output.sequences,
        local = rules.filter_length.output.sequences
    output:
        sequences =  os.path.join(
            build_dir,
            'data',
            '{subtype}_{segment}_concatenated.fasta')
    shell: 
        """
        cat {input.reference} {input.local} > {output.sequences}
        """

"""This rule merges metadata for local and reference samples."""
rule merge_metadata:
    message:
        """
        Merging the metadata for the reference and the local files.\n
        Reference FASTA: {input.reference}\n
        Local FASTA: {input.local}
        """
    input:
        reference = rules.load_lineages.output.metadata,
        local = rules.filter_length.output.metadata
    output:
        metadata = os.path.join(
            build_dir,
            'data',
            '{subtype}_{segment}_merged.tsv')
    shell:
        """
        augur merge \
            --metadata \
                REFERENCE={input.reference} \
                LOCAL={input.local} \
            --output-metadata {output.metadata}
        """

## ALIGN
rule align_lineage:
    message: "Aligning lineage reference sequences with local sequences"
    input:
        fasta = rules.cat_fastas.output.sequences,
        reference_file = lambda wc: ref_files[wc.subtype][wc.segment]['FASTA']
    output:
        alignment = os.path.join(
            build_dir,
            'align',
            '{subtype}_{segment}_alignment.fasta')
    shell:
        """
        augur align \
            --sequences {input.fasta} \
            --reference-sequence {input.reference_file} \
            --output {output.alignment} \
            --nthreads 4 \
            --remove-reference
        """

## BUILD TREE

### RAW TREE
"""This rule builds the initial tree."""
rule tree_lineage:
    message: "Building lineage tree"
    input:
        alignment = rules.align_lineage.output.alignment
    output:
        tree = os.path.join(
            build_dir,
            'tree',
            '{subtype}_{segment}_raw.nwk')
    params:
        method = "iqtree"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method {params.method} \
            --nthreads 4
        """

### Refine
"""This rule refines the tree. In this case, we basically
just use it to assign internal node names, since we aren't
going to assign time to the tree."""
rule refine_lineage:
    message:
        """
        Refining tree and assigning internal node names.
        """
    input:
        tree = rules.tree_lineage.output.tree,
        alignment = rules.align_lineage.output.alignment,
        metadata = rules.merge_metadata.output.metadata
    output:
        tree = os.path.join(
            build_dir,
            'tree',
            '{subtype}_{segment}_refined.nwk'),
        node_data = os.path.join(
            build_dir,
            'tree',
            '{subtype}_{segment}_branch_lengths.json')
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --divergence-units 'mutations-per-site'
        """

### NODES AND TRAITS
"""This rule reconstructs the ancestral node sequences.
It then annotates nucleotide mutations at each branch."""
rule ancestral_lineage:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine_lineage.output.tree,
        alignment = rules.align_lineage.output.alignment,
        reference = lambda wc: ref_files[wc.subtype][wc.segment]['FASTA']
    output:
        nt_muts = os.path.join(
            build_dir,
            'tree',
            '{subtype}_{segment}_nt_muts.json'),
        sequences = os.path.join(
            build_dir,
            'tree',
            '{subtype}_{segment}_ancestral_sequences.fasta')
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.nt_muts} \
            --output-sequences {output.sequences} \
            --root-sequence {input.reference}
        """

"""This rule simply translates the sequence of each gene at each node,
including inferred ancestral nodes."""
rule translate_lineage:
    message: "Translating amino acid sequences and identifying mutations"
    input:
        tree = rules.refine_lineage.output.tree,
        ancestral_json = rules.ancestral_lineage.output.nt_muts,
        reference = lambda wc: ref_files[wc.subtype][wc.segment]['gbk']
    output:
        aa_muts = os.path.join(
            build_dir,
            'tree',
            '{subtype}_{segment}_aa_muts.json'),
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.ancestral_json} \
            --reference-sequence {input.reference} \
            --output-node-data {output.aa_muts}
        """

"""This rule labels HA clades based on defining mutations."""
rule label_clades:
    message: "Labeling clades based on mutations"
    input:
        tree = rules.refine_lineage.output.tree,
        mutations = rules.translate_lineage.output.aa_muts,
        clades = lambda wc: ref_files[wc.subtype][wc.segment]['clade_muts']
    output:
        clade_labels = os.path.join(
            build_dir,
            'tree',
            '{subtype}_{segment}_clade_labels.json'),
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.mutations} \
            --clades {input.clades} \
            --output-node-data {output.clade_labels}
        """

## VISUALIZATION
"""This rule exports the results of the pipeline into JSON
for visualization in auspice."""
rule export_lineage:
    message: "Exporting lineage JSON files for for auspice"
    input:
        tree = rules.refine_lineage.output.tree,
        metadata = rules.merge_metadata.output.metadata,
        node_data = [rules.refine_lineage.output.node_data,
                     rules.ancestral_lineage.output.nt_muts,
                     rules.translate_lineage.output.aa_muts],
                     #rules.label_clades.output.clade_labels],
        auspice_config = 'input/config/build-lineages/auspice_config.json',
    output:
        auspice_json = os.path.join(
            'output',
            'auspice',
            'lineages_{subtype}_{segment}.json'),
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json} \
            --include-root-sequence-inline
        """