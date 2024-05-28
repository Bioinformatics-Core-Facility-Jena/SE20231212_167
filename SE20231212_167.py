#conda create -n SE20231212_167
#conda activate SE20231212_167
#conda install ipython
#conda install fastqc
#conda install -c bioconda multiqc
#conda install fastp
#conda install hisat2
#conda install -c bioconda gffread
#conda install -c bioconda subread
#conda install -c bioconda bioconductor-deseq2
#conda install seaborn
#conda install pandas
#conda install scikit-learn
#conda install upsetplot
#conda install -c conda-forge r-gprofiler2
#conda install -c bioconda r-rlist
#conda install -c bioconda bioconductor-enrichplot
#conda install -c bioconda bioconductor-dose


import os
from glob import glob

workingDirectory = '' #CHANGE THIS VARIABLE TO THE DESIRED WORKING DIRECTORY

##################################################
##### Quality Control with fastQC and fastp
##################################################

os.makedirs(f'{workingDirectory}/trimmed_data/')

for idx,readFile in enumerate(glob(f'{workingDirectory}/raw_data/*gz')):
    os.system(f'fastqc -t 40 --noextract {readFile}')

os.system(f'multiqc --outdir {workingDirectory}/raw_data/ {workingDirectory}/raw_data/')

for readfile in glob(f'{workingDirectory}/raw_data/*_1.fastq.gz'):
    os.system(f'fastp --in1 {readfile}'\
                    f' --in2 {readfile.replace("_1.fastq.gz", "_2.fastq.gz")}'\
                    f' --out1 {workingDirectory}/trimmed_data/{os.path.basename(readfile).replace(".fastq.gz", ".trimmed.fq.gz")}'\
                    f' --out2 {workingDirectory}/trimmed_data/{os.path.basename(readfile).replace("_1.fastq.gz", "_2.trimmed.fq.gz")}'\
                    f' --detect_adapter_for_pe --cut_front --cut_tail --n_base_limit=1 --trim_poly_g --trim_poly_x --length_required=15 --average_qual=28 --thread=16 --html={workingDirectory}/trimmed_data/{os.path.basename(readfile).split("_R1_001.fastq.gz")[0]}.html')

for idx,readFile in enumerate(glob(f'{workingDirectory}/trimmed_data/*gz')):
    os.system(f'fastqc -t 40 --noextract {readFile}')

os.system(f'multiqc --outdir {workingDirectory}/trimmed_data/ {workingDirectory}/trimmed_data/')

##################################################
##### Read Mapping onto A. thaliana genome TAIR10.1tair with hisat2
##################################################

### Download A. thaliana genome TAIR10.1
os.makedirs(f'{workingDirectory}/genome/')
os.system(f'curl -o {workingDirectory}/genome/tair10.1.zip -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001735.4/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,SEQUENCE_REPORT&filename=GCF_000001735.4.zip" -H "Accept: application/zip"')
os.system(f'unzip -d {workingDirectory}/genome/ {workingDirectory}/genome/tair10.1.zip')
os.system(f'mv {workingDirectory}/genome/ncbi_dataset/data/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna {workingDirectory}/genome/tair10.1.fna')
os.system(f'mv {workingDirectory}/genome/ncbi_dataset/data/GCF_000001735.4/genomic.gff {workingDirectory}/genome/tair10.1.gff')

### convert gff to gtf with gffread
os.system(f'gffread --keep-genes -T {workingDirectory}/genome/tair10.1.gff -o {workingDirectory}/genome/tair10.1.gtf')

# the gene_id entry has to be the first entry in the 9th column of each line in the gtf file
cache = open(f'{workingDirectory}/genome/tair10.1.gtf').readlines()
with open(f'{workingDirectory}/genome/tair10.1.gtf', 'w') as outFile:
    for line in cache:
        if 'gene_id' in line:
            currentGeneID = line.split('\t')[8].split('gene_id')[1].split('-')[1].split('"')[0]
            newLine = line.split('transcript_id "')[0] + f'gene_id "{currentGeneID}"; transcript_id' + line.split('transcript_id')[1]
            outFile.write(newLine)

### Indexing A. thaliana genome TAIR10.1 with hisat2
os.system(f'hisat2-build -p 40 {workingDirectory}/genome/tair10.1.fna {workingDirectory}/genome/tair10.1')

### Read Mapping onto A. thaliana genome TAIR10.1 with hisat2
os.makedirs(f'{workingDirectory}/mapping/')
for idx,readFile in enumerate(glob(f'{workingDirectory}/trimmed_data/*_1.trimmed.fq.gz')):
    os.system(f'hisat2 --new-summary --summary-file {workingDirectory}/mapping/{os.path.basename(readFile).replace("_1.trimmed.fq.gz", ".summary")} -p 60 -x {workingDirectory}/genome/tair10.1 -1 {readFile} -2 {readFile.replace("_1.trimmed.", "_2.trimmed.")} -S {workingDirectory}/mapping/{os.path.basename(readFile).replace("_1.trimmed.fq.gz", ".sam")}')
    ### Convert SAM to sorted BAM and index
    samFile = f'{workingDirectory}/mapping/{os.path.basename(readFile).replace("_1.trimmed.fq.gz", ".sam")}'
    os.system(f'samtools view -@ 60 -b {samFile} | samtools sort -@ 40 -o {samFile.replace(".sam", ".bam")}')
    os.system(f'samtools index {samFile.replace(".sam", ".bam")}')
    ### Remove SAM files
    os.system(f'rm {samFile}')

##################################################
##### Read Counting with featureCounts
##################################################

os.makedirs(f'{workingDirectory}/counts/')
bamFiles = sorted(glob(f'{workingDirectory}/mapping/*.bam'))
os.system(f'featureCounts -p --countReadPairs -O -M -T 40 -a {workingDirectory}/genome/tair10.1.gtf -o {workingDirectory}/counts/all_samples_counts.txt {" ".join(bamFiles)} > {workingDirectory}/counts/read_assignment_stats.txt 2>&1')

##################################################
##### Perform differential expression analysis with DESeq2
##################################################

os.makedirs(f'{workingDirectory}/deseq2/')
os.system(f'Rscript {workingDirectory}/deseq2_comparison.r')

##################################################
##### Extend the deseq2 results with gene names from TAIR10.1
##################################################

genesDict = {}
for line in open(f'{workingDirectory}/genome/tair10.1.gff'):
    if line[0] != '#':
        line = line.split('\t')
        if line[2] == 'gene':
            geneID = line[-1].split('ID=')[1].split(';')[0]
            if 'Name=' in line[-1]:
                geneName = line[-1].split('Name=')[1].split(';')[0]
            else:
                geneName = 'NA'
            genesDict[geneID] = geneName

for resultFile in glob(f'{workingDirectory}/deseq2/*/result.csv'):
    with open(resultFile) as inFile:
        with open(resultFile.replace('.csv', '.extended.csv'), 'w') as outFile:
            for idx,line in enumerate(inFile):
                if idx == 0:
                    outFile.write(line.strip() + ',"gene_name"\n')
                else:
                    outFile.write(f'{line.strip()},{genesDict.get(line.split(",")[0][1:-1], "NA")}\n')
