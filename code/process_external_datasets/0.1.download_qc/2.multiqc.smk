"""
Snakemake module for quality control of sequencing reads.

This module handles:
- Running FastQC on individual samples
- Aggregating FastQC reports with MultiQC
- Creating trigger files for dynamic sample discovery
"""

import os
import yaml
import pandas as pd

# Configuration variables
conda_dir = config['conda_dir']
scripts_dir = config['scripts_dir']

wildcard_constraints:
    study_name = '|'.join(config["study_names"])

##############################################
#                FUNCTIONS                   #
##############################################

def get_fastqc(wildcards):
    """
    Generate list of FastQC report paths for all samples in a study.
    
    Args:
        wildcards: Snakemake wildcards containing study_name and strand
        
    Returns:
        list: Paths to FastQC data files for all samples
        
    Note:
        Uses checkpoint to dynamically determine samples after qc_get_samples completes
    """
    with open(checkpoints.qc_get_samples.get(**wildcards).output[0], "r") as file:
        samples = yaml.safe_load(file)[wildcards.study_name]
    
    return expand(
        config['study_dir'] + 'QualityControl/{study_name}/fastqc/{strand}/{sample}/{sample}_{strand}_fastqc/fastqc_data.txt',
        study_name=[wildcards.study_name],
        sample=samples,
        strand=[wildcards.strand]
    )


##############################################
#                  RULES                     #
##############################################

rule multiqc_all:
    """
    Master rule to generate MultiQC reports for all studies and strands.
    
    This rule aggregates QC results across all configured studies and read strands.
    """
    input:
        multiqc_report = expand(
            config['study_dir'] + 'QualityControl/{study_name}/multiqc_{strand}/multiqc_report.html',
            study_name=config['study_names'],
            strand=config['strand']
        )


checkpoint qc_get_samples:
    """
    Create trigger YAML file with list of samples for quality control.
    
    This checkpoint:
        - Reads study metadata
        - Filters samples if needed (AMPLICON strategy for PRJNA383868)
        - Creates a YAML file mapping study_name to sample list
        - Enables dynamic sample discovery for downstream FastQC rules
    """
    input:
        config['study_dir'] + 'raw/{study_name}/{study_name}.csv'
    output:
        config["service_config"] + "multiqc/{study_name}.yaml"
    params:
        study_name = "{study_name}"
    run:
        # Read metadata
        df = pd.read_csv(input[0])
        
        # Special filtering for PRJNA383868 study
        if params.study_name == 'PRJNA383868':
            df = df.loc[df.library_strategy.values == 'AMPLICON', :]
        
        # Extract sample IDs
        samples = df.run_accession.values.tolist()
        
        # Create YAML configuration
        config_dict = {params.study_name: samples}
        
        with open(output[0], "w") as file:
            yaml.dump(config_dict, file, default_flow_style=False)


rule fastqc:
    """
    Run FastQC quality control on a single FASTQ file.
    
    FastQC generates:
        - Quality metrics for sequencing reads
        - Per-base quality scores
        - Sequence content analysis
        - Overrepresented sequences
        - Adapter content detection
    
    The zipped output is automatically extracted for downstream processing.
    """
    input:
        config['study_dir'] + 'raw/{study_name}/raw_reads/{sample}/{sample}_{strand}.fastq.gz'
    output:
        config['study_dir'] + 'QualityControl/{study_name}/fastqc/{strand}/{sample}/{sample}_{strand}_fastqc/fastqc_data.txt'
    params:
        folder = config['study_dir'] + 'QualityControl/{study_name}/fastqc/{strand}/{sample}/',
        zipped = config['study_dir'] + 'QualityControl/{study_name}/fastqc/{strand}/{sample}/{sample}_{strand}_fastqc.zip'
    threads: 10
    # group: lambda wildcards: f"fastqc_{hash(wildcards.sample) % 100}"
    conda:
        f'{conda_dir}/conda_env.qc.yaml'
    shell:
        """
        fastqc -t {threads} -o {params.folder} {input} && \
        unzip -o {params.zipped} -d {params.folder}
        """


rule multiqc:
    """
    Aggregate FastQC reports into a single MultiQC HTML report.
    
    MultiQC:
        - Combines all FastQC results for a study/strand
        - Generates interactive HTML report with plots
        - Provides summary statistics across all samples
        - Highlights samples with quality issues
    
    Note:
        One report is generated per study per strand (forward/reverse reads analyzed separately)
    """
    input:
        get_fastqc
    output:
        config['study_dir'] + 'QualityControl/{study_name}/multiqc_{strand}/multiqc_report.html'
    params:
        fastqc_folder = config['study_dir'] + 'QualityControl/{study_name}/fastqc/{strand}',
        folder = config['study_dir'] + 'QualityControl/{study_name}/multiqc_{strand}'
    conda:
        f'{conda_dir}/conda_env.qc.yaml'
    shell:
        """
        multiqc {params.fastqc_folder} -o {params.folder}
        """