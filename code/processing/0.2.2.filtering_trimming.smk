def get_samples():
    '''
    Get the list of all samples in the dataset
    '''
    samples = os.listdir(config['study_dir']+'reads/raw/')
    return samples

rule all:
    '''
    Request filtering and trimming for all samples in the dataset.
    '''
    input: expand(config['study_dir']+'reads/filtered/{sample}/{sample}_{strand}.fastq.gz', sample = get_samples(), strand = ['R1', 'R2'])

rule dada_filter_and_trim:
    '''
    Run dada2 filter_and_trim function for each sample.
    Creates a new folder - "filtered" in "reads" directory.
    Starts after initial FastQC report is generated.
    '''
    input:
        qc_trigger = config['study_dir']+'QCcontrol/raw/multiqc/multiqc_report_R1.html',
        R1 = config['study_dir']+'reads/raw/{sample}/{sample}_R1.fastq.gz',
        R2 = config['study_dir']+'reads/raw/{sample}/{sample}_R2.fastq.gz'
    output:
        R1 = config['study_dir']+'reads/filtered/{sample}/{sample}_R1.fastq.gz',
        R2 = config['study_dir']+'reads/filtered/{sample}/{sample}_R2.fastq.gz'
    params:
        script_dir = config['script_dir'],
        dada2_params = config['dada2_params']
    conda: config['conda']
    shell: "Rscript {params.script_dir}/dada_filter_trim.r {input} {output} {params.dada2_params}"