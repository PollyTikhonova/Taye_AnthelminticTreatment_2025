##### A DADA2 pipeline for processing filtered reads and creates phyloseq objects.
# This is a general pipeline, that might allows to accomodate single-end sequencing data.
# The pipeline consists of the following steps:
# 1. learn_errors: run a dada2 error model.
# 2. derep: dereplicate reads.
# 3. dada: run the core dada2 algorithm of ASV identification separately for forward and reverse strands.
# 4. merge_pairs: merge dada2 results, recieved from two strands.
# 5. make_seqtab: make an ASV table
# 6. seqtab_nochim: remove bimeras
# 7. assign_taxa: assign taxonomy
# 8. convert_rds: create a phyloseq object and incorporates metadata.
# 9. transform_rds: filter out non-Bacteria or Archaea ASVs, and mitochondria ASVs. 
#     Agglomerates ASV at the genus and family levels.
#     The resulting phyloseq objects are saved in the .rds format.
#     In addition, for compatibility with other programming languages, 
#     the each part of the phyloseq object data (otu_table, tax_table, sample_data) is saved in the form of csv tables.
#     The resulting files are:
#         - nohost_asv: non-agglomerated (ASVs);
#         - nohost_genus: agglomeratd at the genus level;
#         - nohost_family: agglomeratd at the family level.

conda_dir = config["conda_dir"]

def get_seqtab_input(wildcards):
    '''
    A flexibility function, that identifies if the dataset consists 
    if single-end or paired-end sequencing reads.
    In case of single-end reads, it allows the pipeline to skip the merging step and
    directs the results dada rule to the make_seqtab rule.
    '''
    import pandas as pd
    # reads the paths of the samples table
    sample_table = checkpoints.sample_table.get(**wildcards).output
    # reads the samples table
    sample_table = pd.read_csv(str(sample_table), sep='\t')
    # If the samples table has a column for reverse reads, treats the dataset as paired-end.
    # Otherwise, it considers the dataset single-end.
    if 'R2' in sample_table.columns.values.tolist():
        return config['study_dir']+'dada2/{study_name}/mergers.rds'
    else:
        return config['study_dir']+'dada2/{study_name}/dada_R1.rds'

rule all:
    '''
    Requests the unfiltered and filtered phyloseq objects as the result of the pipeline.
    Filtered objects are requested with additional agglomerations at genus and family levels.
    Filtered object are generated in phyloseq.rds objects, as well as .csv tables.
    '''
    input:         
        phyloseq_raw = expand(config['study_dir']+"phyloseq/{study_name}/phyloseq.rds", study_name=config['study_names']), 
        phyloseq_nohost = expand(config['study_dir']+"phyloseq/{study_name}/nohost_{agglom}/phyloseq.rds", study_name=config['study_names'], agglom = ['asv', 'genus', 'family'])
   
rule learn_errors:
    '''
    Runs a error model for all samples. In case of paired end sequencing, 
    each strand is processed separately.
    Creates one object per strand for all samples in the dataset.
    '''
    input: config["study_dir"] + "processed/trimmed/{study_name}.csv"
    output: config['study_dir']+'dada2/{study_name}/err_{strand}.rds'
    log: config['study_dir']+'dada2/{study_name}/logs/err_{strand}.log'
    conda:
        f"{conda_dir}/conda_env.dada2.yaml"
    params:
        strand = lambda w: w.strand,
        script_dir = config['script_dir']
    shell: "Rscript {params.script_dir}/learn_err.r {input} {params.strand} {output} > {log}"
    
rule derep:
    '''
    Dereplicates the reads for each sample in the dataset.
    Creates one object per strand for all samples in the dataset.
    '''
    input: config['study_dir']+'dada2/{study_name}/sample_table.csv',
    output: config['study_dir']+'dada2/{study_name}/derep_{strand}.rds'
    conda:
        f"{conda_dir}/conda_env.dada2.yaml"
    params:
        script_dir = config['script_dir'],
        strand = lambda w: w.strand
    shell: "Rscript {params.script_dir}/derep.r {input} {params.strand} {output}"
    
rule dada:
    '''
    dada2 main algorithm: utilizes a error model to identify ASVs on dereplicated data.
    Creates one object per strand for all samples in the dataset.
    '''
    input: 
        err_file = config['study_dir']+'dada2/{study_name}/err_{strand}.rds',
        derep_dependency = config['study_dir']+'dada2/{study_name}/derep_{strand}.rds'
    output: config['study_dir']+'dada2/{study_name}/dada_{strand}.rds'
    conda:
        f"{conda_dir}/conda_env.dada2.yaml"
    params:
        script_dir = config['script_dir'],
        pool = config['pool']
    shell: "Rscript {params.script_dir}/dada.r {input.derep_dependency} {input.err_file} {output} {params.pool}"
    
rule merge_pairs:
    '''
    Merges dada ASVs from two strands into one object.
    '''
    input: 
        dada_file = expand(config['study_dir']+'dada2/{study_name}/dada_{strand}.rds', strand = ['R1', 'R2']),
        derep_file = expand(config['study_dir']+'dada2/{study_name}/derep_{strand}.rds', strand = ['R1', 'R2'])
    output: config['study_dir']+'dada2/{study_name}/mergers.rds'
    conda:
        f"{conda_dir}/conda_env.dada2.yaml"
    params:
        script_dir = config['script_dir']
    shell: "Rscript {params.script_dir}/merge_pairs.r {input.dada_file} {input.derep_file} {output};"

rule make_seqtab:
    '''
    Creates an ASV table.
    '''
    input: get_seqtab_input
    output: config['study_dir']+'dada2/{study_name}/seqtab.rds'
    conda:
        f"{conda_dir}/conda_env.dada2.yaml"
    params:
        script_dir = config['script_dir']
    log: config['study_dir']+'dada2/{study_name}/logs/make_seqtab.log'
    shell: "Rscript {params.script_dir}/make_seqtab.r {input} {output}"

rule seqtab_nochim:
    '''
    Removes bimeras.
    '''
    input: config['study_dir']+'dada2/{study_name}/seqtab.rds'
    output: config['study_dir']+'dada2/{study_name}/seqtab_nochim.rds'
    conda:
        f"{conda_dir}/conda_env.dada2.yaml"
    params:
        script_dir = config['script_dir']
    log: config['study_dir']+'dada2/{study_name}/logs/seqtab_nochim.log'
    shell: "Rscript {params.script_dir}/seqtab_nochim.r {input} {output}"

rule assign_taxa:
    '''
    Assigns taxonomy using SILVA database.
    '''
    input: config['study_dir']+'dada2/{study_name}/seqtab_nochim.rds'
    output: config['study_dir']+'dada2/{study_name}/taxa.rds'
    conda:
        f"{conda_dir}/conda_env.dada2.yaml"
    params:
        script_dir = config['script_dir'],
        db_silva = config['db_silva'],
        db_species = config['db_species']
    log: config['study_dir']+'dada2/{study_name}/logs/assign_taxa.log'
    shell: "Rscript {params.script_dir}/assign_taxa.r {input} {output} {params.db_silva} {params.db_species}"

rule convert_rds:
    '''
    Creates a phyloseq object that contains:
    read counts (otu_table), taxonomy (tax_table), and provided metadata (sample_data).
    '''
    input: 
        seqtab_nochim = config['study_dir']+'dada2/{study_name}/seqtab_nochim.rds',
        taxa = config['study_dir']+'dada2/{study_name}/taxa.rds',
        metadata = config['metadata']
    output: config['study_dir']+"phyloseq/{study_name}/phyloseq.rds"
    params:
        script_dir = config['script_dir'],
    log: config['study_dir']+'logs/dada2/{study_name}/convert.log'
    conda:
        f"{conda_dir}/conda_env.dada2.yaml"
    shell: "Rscript {params.script_dir}/convert.r {input} {output}"

rule transform_rds:
    '''
    Removes ASVs:
     - not mapped to Bacteria or Archaea Kingdoms;
     - mapped to Mitochondria family.
    Agglomerates ASVs at genus and family levels.
    Saves phyloseq object in .rds format,
    as well we separate .csv dataframes for otu_table, sample_data, tax_table.
    '''
    input: config['study_dir']+"phyloseq/{study_name}/phyloseq.rds"
    output: 
        ps = config['study_dir']+"phyloseq/{study_name}/nohost_{agglom}/phyloseq.rds",
        otu = config['study_dir']+"phyloseq/{study_name}/nohost_{agglom}/otu_table.csv",
        taxa = config['study_dir']+"phyloseq/{study_name}/nohost_{agglom}/tax_table.csv",
        sample_data = config['study_dir']+"phyloseq/{study_name}/nohost_{agglom}/sample_data.csv"
    params:
        agglom = "{agglom}",
        script_dir = config['script_dir']
    conda:
        f"{conda_dir}/conda_env.dada2.yaml"
    shell: "Rscript {params.script_dir}/transform.r {input} {output} {params.agglom}"

