rule all:
    '''
    Request a MultiQC report
    '''
    input:
        mutliqc_report = expand(config['study_dir']+'QCcontrol/{reads}/multiqc/multiqc_report_{strand}.html', reads=['raw'], strand=config['strand'])

def get_fastqc(wildcards):
    '''
    Create a list of all sample names and reqiest a FastQC report for each of them
    '''
    samples = os.listdir(config['study_dir']+'reads/raw/')
    return expand(config['study_dir']+'QCcontrol/{{reads}}/fastqc/{{strand}}/{sample}/{sample}_{strand}_fastqc/fastqc_data.txt', sample=samples, strand=config['strand'])

rule fastq:
    '''
    Run FastQC for a sample and unzip the resulting file
    '''
    input: config['study_dir']+'reads/{reads}/{sample}/{sample}_{strand}.fastq.gz'
    output: config['study_dir']+'QCcontrol/{reads}/fastqc/{strand}/{sample}/{sample}_{strand}_fastqc/fastqc_data.txt'
    params: 
        folder = config['study_dir']+'QCcontrol/{reads}/fastqc/{strand}/{sample}/',
        zipped = config['study_dir']+'QCcontrol/{reads}/fastqc/{strand}/{sample}/{sample}_{strand}_fastqc.zip'
    threads: 10
    conda: config['conda']
    shell: "fastqc -t {threads} -o {params.folder} {input}; unzip -o {params.zipped} -d {params.folder}"

rule multiqc:
    '''
    Gather all FastQC reports for the same strand and generate a MultiQC html report.
    '''
    input:  get_fastqc
    output: config['study_dir']+'QCcontrol/{reads}/multiqc/multiqc_report_{strand}.html'
    params:
        fastqc_folder = config['study_dir']+'QCcontrol/{reads}/fastqc/{strand}',
        folder = config['study_dir']+'QCcontrol/{reads}/{strand}',
        file = config['study_dir']+'QCcontrol/{reads}/{strand}/multiqc_report.html'
    conda: config['conda']
    shell: "multiqc {params.fastqc_folder} -o {params.folder}; mv {params.file} {output}"