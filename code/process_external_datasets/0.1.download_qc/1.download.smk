"""
Snakemake module for downloading sequencing reads from SRA.

This module handles:
- Downloading metadata for studies
- Downloading paired-end and single-end reads
- Recording batch information for downstream processing
"""

import pandas as pd
import os

# Rule execution order: single-end downloads take precedence over paired-end
ruleorder: download_reads_single > download_reads_paired

wildcard_constraints:
    study_name = '|'.join(config["study_names"])

# Configuration variables
conda_dir = config['conda_dir']
scripts_dir = config['scripts_dir']

##############################################
#                FUNCTIONS                   #
##############################################


def sample_files(wildcards):
    """
    Generate list of expected FASTQ files based on study metadata.
    
    Args:
        wildcards: Snakemake wildcards containing study_name
        
    Returns:
        list: Paths to expected FASTQ files for all samples in the study
        
    Note:
        - Determines single vs paired-end based on library_layout
    """
    study_table = checkpoints.download_metadata.get(**wildcards).output[0]
    study_name = wildcards.study_name

    df = pd.read_csv(study_table)

    
    samples = df.run_accession.values

    # Determine if reads are paired-end or single-end
    if 'PAIRED' in df.library_layout.values:
        strands = ['2']
    else:
        strands = ['1']
    
    return expand(
        config['study_dir'] + 'raw/{study_name}/raw_reads/{sample}/{sample}_{strand}.fastq.gz',
        study_name=[study_name],
        sample=samples,
        strand=strands
    )

##############################################
#                  RULES                     #
##############################################

rule download_all:
    """
    Master rule to download all metadata and batch files for all studies.
    """
    input:
        metadata = expand(
            config['study_dir'] + 'raw/{study_name}/{study_name}.csv',
            study_name=config['study_names']
        ),
        batch_files = expand(
            config['study_dir'] + 'raw/{study_name}/batches.csv',
            study_name=config['study_names']
        )


checkpoint download_metadata:
    """
    Download metadata for a study from SRA.
    
    This checkpoint allows dynamic determination of samples after metadata is retrieved.
    """
    output:
        config['study_dir'] + 'raw/{study_name}/{study_name}.csv'
    params:
        study_name = '{study_name}',
        study_dir = config['study_dir'] + 'raw/{study_name}/'
    conda:
        f'{conda_dir}/download_trim_python.yaml'
    shell:
        f"""
        python {scripts_dir}/download_trim/download.py {{params}} {{output}}
        """


rule download_reads_paired:
    """
    Download paired-end sequencing reads from SRA.
    
    Steps:
        1. Prefetch SRA data
        2. Validate downloaded data
        3. Extract FASTQ files (split into forward/reverse)
        4. Clean up intermediate SRA files
    """
    input:
        ancient(config['study_dir'] + 'raw/{study_name}/{study_name}.csv')
    output:
        forward_ = config['study_dir'] + 'raw/{study_name}/raw_reads/{sample}/{sample}_1.fastq.gz',
        reverse_ = config['study_dir'] + 'raw/{study_name}/raw_reads/{sample}/{sample}_2.fastq.gz'
    params:
        folder = config['study_dir'] + 'raw/{study_name}/raw_reads/{sample}/',
        sra_id = '{sample}'
    # group:
    #     lambda wildcards: f"download_reads_{hash(wildcards.sample) % 100}"
    conda:
        f'{conda_dir}/download_trim.yaml'
    shell:
        """
        prefetch --output-directory {params.folder} {params.sra_id} && \
        vdb-validate {params.folder}{params.sra_id} && \
        fastq-dump --split-files --gzip --outdir {params.folder} {params.folder}{params.sra_id} && \
        rm -r {params.folder}{params.sra_id}
        """


rule download_reads_single:
    """
    Download single-end sequencing reads from SRA.
    
    Steps:
        1. Prefetch SRA data
        2. Validate downloaded data
        3. Extract FASTQ file
        4. Rename to _1.fastq.gz for consistency
        5. Clean up intermediate SRA files
    """
    input:
        ancient(config['study_dir'] + 'raw/{study_name}/{study_name}.csv')
    output:
        config['study_dir'] + 'raw/{study_name}/raw_reads/{sample}/{sample}_1.fastq.gz'
    params:
        folder = config['study_dir'] + 'raw/{study_name}/raw_reads/{sample}/',
        sra_id = '{sample}'
    # group:
    #     lambda wildcards: f"download_reads_{hash(wildcards.sample) % 100}"
    conda:
        f'{conda_dir}/download_trim.yaml'
    shell:
        """
        prefetch --output-directory {params.folder} {params.sra_id} && \
        vdb-validate {params.folder}{params.sra_id} && \
        fastq-dump --gzip --outdir {params.folder} {params.folder}{params.sra_id} && \
        mv {params.folder}{params.sra_id}.fastq.gz {params.folder}{params.sra_id}_1.fastq.gz && \
        rm -r {params.folder}{params.sra_id}
        """


rule record_batches:
    """
    Create batch information file for downloaded samples.
    
    This rule:
        - Filters metadata to only include successfully downloaded samples
        - Extracts batch information (using 'sequencing run' column for PRJEB13747)
        - Records strand information (single vs paired-end)
        - Outputs a CSV file for downstream batch effect correction
    """
    input:
        metadata = config['study_dir'] + 'raw/{study_name}/{study_name}.csv',
        downloaded_files = sample_files
    output:
        config['study_dir'] + 'raw/{study_name}/batches.csv'
    params:
        study_name = '{study_name}',
        files_directory = config['study_dir'] + 'raw/{study_name}/raw_reads/'
    run:
        # Read metadata
        data = pd.read_csv(input.metadata)
        data.index = data.run_accession.values
        
        # Get list of successfully downloaded samples
        samples = os.listdir(params.files_directory)
        data = data.loc[samples, :]
        data.index = list(range(len(data)))

        # Determine strand information
        if 'PAIRED' in data.library_layout.values:
            strands = ['1;2']
        else:
            strands = ['1']

        # Create batch dataframe
        # Default: all samples in single batch
        df_batches = pd.DataFrame({
            'run_accession': data['run_accession'].values,
            'batch': ['1'] * len(data),
            'study_accession': data['study_accession'].values,
            'strands': strands * len(data)
        })
        
        # Save batch information
        df_batches.to_csv(output[0], index=False)