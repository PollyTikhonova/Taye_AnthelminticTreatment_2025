##### A DADA2 pipeline for processing filtered reads and creates phyloseq objects.
# This is a general pipeline, that might allows to accomodate single-end sequencing data.
# The pipeline consists of the following steps:
# 1. sample_table: create a .csv table that contains a list of all samples in the dataset.
# 2. learn_errors: run a dada2 error model.
# 3. derep: dereplicate reads.
# 4. dada: run the core dada2 algorithm of ASV identification separately for forward and reverse strands.
# 5. merge_pairs: merge dada2 results, recieved from two strands.
# 6. make_seqtab: make an ASV table
# 7. seqtab_nochim: remove bimeras
# 8. assign_taxa: assign taxonomy
# 9. convert_rds: create a phyloseq object and incorporates metadata.
# 10. transform_rds: filter out non-Bacteria or Archaea ASVs, and mitochondria ASVs. 
#     Agglomerates ASV at the genus and family levels.
#     The resulting phyloseq objects are saved in the .rds format.
#     In addition, for compatibility with other programming languages, 
#     the each part of the phyloseq object data (otu_table, tax_table, sample_data) is saved in the form of csv tables.
#     The resulting files are:
#         - nohost_asv: non-agglomerated (ASVs);
#         - nohost_genus: agglomeratd at the genus level;
#         - nohost_family: agglomeratd at the family level.

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
        return config['study_dir']+'dada2/mergers.rds'
    else:
        return config['study_dir']+'dada2/dada_R1.rds'

rule all:
    '''
    Requests the unfiltered and filtered phyloseq objects as the result of the pipeline.
    Filtered objects are requested with additional agglomerations at genus and family levels.
    Filtered object are generated in phyloseq.rds objects, as well as .csv tables.
    '''
    input:         
        phyloseq_raw = config['study_dir']+"phyloseq/phyloseq.rds",
        phyloseq_nohost = expand(config['study_dir']+"phyloseq/nohost_{agglom}/phyloseq.rds", agglom = ['asv', 'genus', 'family'])

checkpoint sample_table:
    '''
    Creates a sample table that has two columns in case of single-end data: sample, R1;
    and three columns in case of paired-end data: sample, R1, R2.
    "sample" column contains the names of the samples in the dataset, 
    "R1" and "R2" columns contain paths to the correponding sample files, if existent.
    '''
    input: config['study_dir']+'QCcontrol/filtered/multiqc/multiqc_report_R1.html',
    output: config['study_dir']+'dada2/sample_table.csv'
    log: config['study_dir']+'dada2/logs/sample_table.log'
    conda: config['conda']
    params:
        script_dir = config['script_dir'],
        folder = config['study_dir']+'reads/filtered/'
    shell: "Rscript {params.script_dir}/create_sample_table.r {params.folder} {output}"
    
rule learn_errors:
    '''
    Runs a error model for all samples. In case of paired end sequencing, 
    each strand is processed separately.
    Creates one object per strand for all samples in the dataset.
    '''
    input: config['study_dir']+'dada2/sample_table.csv'
    output: config['study_dir']+'dada2/err_{strand}.rds'
    log: config['study_dir']+'dada2/logs/err_{strand}.log'
    conda: config['conda']
    params:
        strand = lambda w: w.strand,
        script_dir = config['script_dir']
    shell: "Rscript {params.script_dir}/learn_err.r {input} {params.strand} {output} > {log}"
    
rule derep:
    '''
    Dereplicates the reads for each sample in the dataset.
    Creates one object per strand for all samples in the dataset.
    '''
    input: config['study_dir']+'dada2/sample_table.csv',
    output: config['study_dir']+'dada2/derep_{strand}.rds'
    conda: config['conda']
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
        err_file = config['study_dir']+'dada2/err_{strand}.rds',
        derep_dependency = config['study_dir']+'dada2/derep_{strand}.rds'
    output: config['study_dir']+'dada2/dada_{strand}.rds'
    conda: config['conda']
    params:
        script_dir = config['script_dir'],
        pool = config['pool']
    shell: "Rscript {params.script_dir}/dada.r {input.derep_dependency} {input.err_file} {output} {params.pool}"
    
rule merge_pairs:
    '''
    Merges dada ASVs from two strands into one object.
    '''
    input: 
        dada_file = expand(config['study_dir']+'dada2/dada_{strand}.rds', strand = ['R1', 'R2']),
        derep_file = expand(config['study_dir']+'dada2/derep_{strand}.rds', strand = ['R1', 'R2'])
    output: config['study_dir']+'dada2/mergers.rds'
    conda: config['conda']
    params:
        script_dir = config['script_dir']
    shell: "Rscript {params.script_dir}/merge_pairs.r {input.dada_file} {input.derep_file} {output};"

rule make_seqtab:
    '''
    Creates an ASV table.
    '''
    input: get_seqtab_input
    output: config['study_dir']+'dada2/seqtab.rds'
    conda: config['conda']
    params:
        script_dir = config['script_dir']
    log: config['study_dir']+'dada2/logs/make_seqtab.log'
    shell: "Rscript {params.script_dir}/make_seqtab.r {input} {output}"

rule seqtab_nochim:
    '''
    Removes bimeras.
    '''
    input: config['study_dir']+'dada2/seqtab.rds'
    output: config['study_dir']+'dada2/seqtab_nochim.rds'
    conda: config['conda']
    params:
        script_dir = config['script_dir']
    log: config['study_dir']+'dada2/logs/seqtab_nochim.log'
    shell: "Rscript {params.script_dir}/seqtab_nochim.r {input} {output}"

rule assign_taxa:
    '''
    Assigns taxonomy using SILVA database.
    '''
    input: config['study_dir']+'dada2/seqtab_nochim.rds'
    output: config['study_dir']+'dada2/taxa.rds'
    conda: config['conda']
    params:
        script_dir = config['script_dir'],
        db_silva = config['db_silva'],
        db_species = config['db_species']
    log: config['study_dir']+'dada2/logs/assign_taxa.log'
    shell: "Rscript {params.script_dir}/assign_taxa.r {input} {output} {params.db_silva} {params.db_species}"

rule convert_rds:
    '''
    Creates a phyloseq object that contains:
    read counts (otu_table), taxonomy (tax_table), and provided metadata (sample_data).
    '''
    input: 
        seqtab_nochim = config['study_dir']+'dada2/seqtab_nochim.rds',
        taxa = config['study_dir']+'dada2/taxa.rds',
        metadata = config['metadata']
    output: config['study_dir']+"phyloseq/phyloseq.rds"
    params:
        script_dir = config['script_dir'],
    log: config['study_dir']+'logs/dada2/convert.log'
    conda: config['conda']
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
    input: config['study_dir']+"phyloseq/phyloseq.rds"
    output: 
        ps = config['study_dir']+"phyloseq/nohost_{agglom}/phyloseq.rds",
        otu = config['study_dir']+"phyloseq/nohost_{agglom}/otu_table.csv",
        taxa = config['study_dir']+"phyloseq/nohost_{agglom}/tax_table.csv",
        sample_data = config['study_dir']+"phyloseq/nohost_{agglom}/sample_data.csv"
    params:
        agglom = "{agglom}",
        script_dir = config['script_dir']
    conda: config['conda']
    shell: "Rscript {params.script_dir}/transform.r {input} {output} {params.agglom}"

