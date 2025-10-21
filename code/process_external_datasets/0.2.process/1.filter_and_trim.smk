"""
Snakemake module for filtering and trimming sequencing reads with DADA2.

This module handles:
- Organizing samples by batch
- Creating per-batch parameter files
- Running DADA2 filterAndTrim on single-end and paired-end reads
- Generating sample tables with trimming statistics
"""

import os
import yaml
import pandas as pd
import numpy as np

# Configuration variables
conda_dir = config["conda_dir"]
scripts_dir = config["scripts_dir"]

# Rule execution order: single-end processing takes precedence over paired-end
ruleorder: dada_filter_and_trim_single > dada_filter_and_trim_paired

wildcard_constraints:
    study_name = '|'.join(config["study_names"])


##############################################
#                FUNCTIONS                   #
##############################################

def get_batchfiles(wildcards):
    """
    Generate list of batch parameter files for a study.
    
    Args:
        wildcards: Snakemake wildcards containing study_name
        
    Returns:
        list: Paths to parameter YAML files for each batch
        
    Note:
        Uses checkpoint to dynamically determine batches after batch info is retrieved
    """
    config_batches = checkpoints.filter_and_trim_batches.get(**wildcards).output[0]
    with open(config_batches, "r") as file:
        batches = list(yaml.safe_load(file)['batches'].keys())
    
    return expand(
        config["params_dir"] + "filter_and_trim/{study_name}_{batch}.yaml",
        study_name=[wildcards.study_name],
        batch=batches
    )

def get_trimming_logs(wildcards):
    """
    Generate list of trimming log files for all samples in a study.
    
    Args:
        wildcards: Snakemake wildcards containing study_name
        
    Returns:
        list: Paths to log files for all samples (using last strand only)
        
    Note:
        For paired-end data, only requests the _2.log file (which implies _1 completed)
    """
    config_batches = checkpoints.filter_and_trim_batches.get(**wildcards).output[0]
    with open(config_batches, "r") as file:
        config_batches_data = yaml.safe_load(file)
    
    # Get the last strand (2 for paired-end, 1 for single-end)
    strands = [config_batches_data['strands'].split(";")[-1]]
    
    # Collect all sample names across batches
    sample_names = []
    for batch, samples in config_batches_data['batches'].items():
        sample_names.extend(samples)
    
    return expand(
        config["study_dir"] + "processed/trimmed/{study_name}/{sample}/{sample}_{strand}.log",
        study_name=[wildcards.study_name],
        sample=sample_names,
        strand=strands
    )


##############################################
#                  RULES                     #
##############################################

rule dada_filter_and_trim_all:
    """
    Master rule to filter and trim all samples across all studies.
    
    Ensures all studies have completed trimming and generated sample tables.
    """
    input:
        config_expansion = expand(
            config["study_dir"] + "processed/trimmed/{study_name}.csv",
            study_name=config["study_names"]
        )


checkpoint filter_and_trim_batches:
    """
    Create YAML configuration mapping samples to batches and strand information.
    
    This checkpoint:
        - Reads batch information from batches.csv
        - Extracts strand configuration (single vs paired-end)
        - Groups samples by batch
        - Creates trigger file for dynamic batch processing
        
    Output structure:
        strands: "1" or "1;2"
        batches:
            batch1: [sample1, sample2, ...]
            batch2: [sample3, sample4, ...]
    """
    input:
        config["study_dir"] + "raw/{study_name}/batches.csv"
    output:
        config["service_config"] + "filter_and_trim/{study_name}.yaml"
    params:
        study_name = "{study_name}"
    run:
        # Read batch information
        data = pd.read_csv(input[0])
        
        # Extract strand information (same for all samples in study)
        strands = str(data.strands.values[0])
        
        # Create batch dictionary
        batch_dict = {'strands': strands, 'batches': {}}
        
        for batch, data_batch in data.groupby('batch'):
            batch_dict['batches'][str(batch)] = data_batch.run_accession.values.tolist()
        
        # Write configuration file
        with open(output[0], "w") as file:
            yaml.dump(batch_dict, file, default_flow_style=False)


rule dada_filter_and_trim_paired:
    """
    Filter and trim paired-end reads using DADA2.
    
    DADA2 filterAndTrim:
        - Removes low-quality bases
        - Trims adapters and primers
        - Filters reads by quality score
        - Removes reads below minimum length
        - Outputs filtered FASTQ files and statistics
        
    Note:
        Uses batch-specific parameters via symlinked parameter files
    """
    input:
        R1 = ancient(config["study_dir"] + "raw/{study_name}/raw_reads/{sample}/{sample}_1.fastq.gz"),
        R2 = ancient(config["study_dir"] + "raw/{study_name}/raw_reads/{sample}/{sample}_2.fastq.gz"),
        dada2_params = config["params_dir"] + "filter_and_trim/{study_name}.yaml"
    output:
        config["study_dir"] + "processed/trimmed/{study_name}/{sample}/{sample}_2.log"
    params:
        R1 = config["study_dir"] + "processed/trimmed/{study_name}/{sample}/{sample}_1.fastq.gz",
        R2 = config["study_dir"] + "processed/trimmed/{study_name}/{sample}/{sample}_2.fastq.gz"
    # group: 
    #     lambda wildcards: f"filter_{hash(wildcards.sample) % 100}"
    conda:
        f"{conda_dir}/dada2.yaml"
    shell:
        f"""
        Rscript {scripts_dir}/download_trim/dada2_filter_trim.r {{input}} {{params}} {{output}}
        """


rule dada_filter_and_trim_single:
    """
    Filter and trim single-end reads using DADA2.
    
    DADA2 filterAndTrim:
        - Removes low-quality bases
        - Trims adapters and primers
        - Filters reads by quality score
        - Removes reads below minimum length
        - Outputs filtered FASTQ files and statistics
        
    Note:
        Uses batch-specific parameters via symlinked parameter files
    """
    input:
        R1 = ancient(config["study_dir"] + "raw/{study_name}/raw_reads/{sample}/{sample}_1.fastq.gz"),
        dada2_params = ancient(get_paramfile)
    output:
        config["study_dir"] + "processed/trimmed/{study_name}/{sample}/{sample}_1.log"
    params:
        R1 = config["study_dir"] + "processed/trimmed/{study_name}/{sample}/{sample}_1.fastq.gz"
    # group: 
    #     lambda wildcards: f"filter_{hash(wildcards.sample) % 100}"
    conda:
        f"{conda_dir}/dada2.yaml"
    shell:
        f"""
        Rscript {scripts_dir}/download_trim/dada2_filter_trim.r {{input}} {{params}} {{output}}
        """


rule create_sample_table:
    """
    Create summary table of trimming statistics for all samples.
    
    This rule:
        - Collects all trimming log files
        - Parses statistics (reads in, reads out, reads filtered)
        - Generates CSV table for quality assessment
        - Used for tracking data loss through filtering pipeline
    """
    input:
        trimming_logs = get_trimming_logs
    output:
        config["study_dir"] + "processed/trimmed/{study_name}.csv"
    params:
        reads_dir = config["study_dir"] + "processed/trimmed/{study_name}/"
    conda:
        f"{conda_dir}/dada2.yaml"
    shell:
        f"""
        Rscript {scripts_dir}/download_trim/create_sample_table.r {{params}} {{output}}
        """