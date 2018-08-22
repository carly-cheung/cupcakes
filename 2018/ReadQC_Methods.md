# Methods to Quality Control Sequencing Reads

## Launch EC2 Instance, Connect to Instance, Create New Directory
```
aegea launch --iam-role S3fromEC2 --ami-tags Name=czbiohub-specops -t m5.4xlarge --duration-hours 2 kalani-cc4

aegea ssh ubuntu@kalani-cc4

cd /mnt/data

mkdir readqc

cd readqc
```
- Info about `specops` image here:https://github.com/czbiohub/packer-images
  - We are using this because FastQC is built into it
- To keep in mind: m5.4xlarge has 16 threads we can use
- We are changing directories into `/mnt/data` because it is a there is 1TB of memory attached to this

## [Skip Step - For Reference] Preprocessing of original reads
```
#Download both read files in one step from seqbot bucket
aws s3 cp --recursive --exclude "*" --include "CMS_089a*" s3://czbiohub-seqbot/fastqs/180523_NS500126_0798_AHWTLTAFXX/ .

#gunzip
gunzip *.gz

#look at read counts in each file
ubuntu@kalani-fastqc4:/mnt/data/fastqc$ awk '{s++}END{print s/4}' CMS_089a_MagBeads_RNA_A_S7_R1_001.fastq
20858066
ubuntu@kalani-fastqc4:/mnt/data/fastqc$ awk '{s++}END{print s/4}' CMS_089a_MagBeads_RNA_A_S7_R2_001.fastq
20858066

#take a subset (5,000,000 reads) of each file and create new file
head -20000000 CMS_089a_MagBeads_RNA_A_S7_R1_001.fastq > CMS_089a_R1_20e6.fastq
head -20000000 CMS_089a_MagBeads_RNA_A_S7_R2_001.fastq > CMS_089a_R2_20e6.fastq

#look at read counts in each file
ubuntu@kalani-fastqc4:/mnt/data/fastqc$ awk '{s++}END{print s/4}' CMS_089a_R1_20e6.fastq
5000000
ubuntu@kalani-fastqc4:/mnt/data/fastqc$ awk '{s++}END{print s/4}' CMS_089a_R2_20e6.fastq
5000000

#Save read subset to S3 bucket for easier access
ubuntu@kalani-fastqc4:/mnt/data/fastqc$ aws s3 cp CMS_089a_R2_20e6.fastq s3://kalani-bucket
```

## Download Data from S3 Bucket
```
aws s3 cp --recursive --exclude "*" --include "CMS_089a*" s3://kalani-bucket .
```

## Find the Trimmomatic Adapter Files

Use the unix `find` command to locate the file path of the adapters called `TruSeq3-PE.fa`:

```
ubuntu@olgabot-cupcakes-kalani:/mnt/data/readqc$ find $HOME -name TruSeq3-PE.fa
/home/ubuntu/anaconda/pkgs/trimmomatic-0.38-0/share/trimmomatic-0.38-0/adapters/TruSeq3-PE.fa
/home/ubuntu/anaconda/share/trimmomatic-0.38-0/adapters/TruSeq3-PE.fa
```

## PriceSeq Filter the starting data
To look at the initial quality of the reads we are using.

PriceSeqFilter Manual: http://derisilab.ucsf.edu/software/price/PriceDocumentation140408/independentQualityFilter.html
```
PriceSeqFilter -fp CMS_089a_R1_20e6.fastq CMS_089a_R2_20e6.fastq -op CMS_089a_R1_20e6_PF.fastq CMS_089a_R2_20e6_PF.fastq -a 12 -rnf 90 -log c -rqf 85 0.98
```
These are the IDSeq parameters:
- a 12: # of threads to use
- rnf 90: 90% of the reads must be called. filters our read pairs if either has an unacceptably high number of uncalled nucleotides (N)
- log c: determines the type of standard output. c = concise stout (default)
- rqf 85 0.98: filters out squences with an unaccepbably high number of low-quality nucleotdies, as well as defined by the provided quality scores (only applies to files whose formats include quality score information)
  - 85: the percent of nucleotides in a read that must be high-quality (determined by the .98)
  - 0.98: the minimium allowed probablility of a nucletide being correct (must be between 0 and 1)  -- determines quality -- Phred scores of 19.91 and above for a base


## Trim reads using Trimmomatic
```
trimmomatic PE -threads 4 -phred33 CMS_089a_R1_20e6.fastq CMS_089a_R2_20e6.fastq CMS_089a_R1_20e6_trim.fastq CMS_089a_R1_20e6_trim_un.fastq CMS_089a_R2_20e6_trim.fastq CMS_089a_R2_20e6_trim_un.fastq ILLUMINACLIP:/home/ubuntu/anaconda/share/trimmomatic-0.38-0/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:30
```
These are the coding carpentry parameters except for min size (20 to 30) adjusted to be compatible with CDHitDup. Additionally, `2:30:10` from manual default.
- PE: paired end
- threads 4: # of threads to use
- phred33: encoding system for quality scores
- ILLUMINACLIP: to clip the Illumina adapters from the input file using the adapter sequences listed in TruSeq3-PE.fa
  - 2:40:15 - Initially Trimmomatic will look for seed matches (16 bases) allowing maximally 2 mismatches. These seeds will be extended and clipped if in the case of paired end reads a score of 30 is reached or in the case of single ended reads a score of 10 (about 17 bases)
- SLIDINGWINDOW: 4:20: to use a sliding window of size 4 that will remove bases if their phred score is below 20
- MINLEN:30: drop an entire read if it is below a specified length

## PriceSeqFilter of Trimmed reads
```
PriceSeqFilter -fp CMS_089a_R1_20e6_trim.fastq CMS_089a_R2_20e6_trim.fastq -op CMS_089a_R1_20e6_trim_PF.fastq CMS_089a_R2_20e6_trim_PF.fastq  -a 12 -rnf 90 -log c -rqf 85 0.98
```

## Take a look at all the files generated in the readqc folder
```
ls -FlhS
```
There should be 10 files and one folder in this.
The -S puts the files in order by size.

## Download all .fastq files onto computer
In your **local terminal**, download all the files generated. Need to perform fastqc on our computer because the program doesn't work on the image itself.

Create a directory in the code folder called fastqc. Or create this directory elsewhere where you know the path.
```
cd ~/code/
mkdir fastqc
```

Download the fastq files generated in instance to your fastqc folder. Make sure you are in the fastqc folder you generated previously. `pwd` to make sure.
```
scp -r <your_ec2_path>:/mnt/data/readqc/ .
```
- You can retrieve your ec2 path by going to your ec2 terminal tab and going to View>Show Inspector. In the 'Running Process' table it should be listed as ssh <your_ec2_path>
  - example of an ec2 path: ` ubuntu@ec2-52-35-162-56.us-west-2.compute.amazonaws.com`
- This will take a while to download so while this happens, open a new tab in terminal that is a local terminal environment.

# Creating a FastQC environment on your computer
## Create an environment and launch it
```
conda update conda
conda create -n fastqc
source activate fastqc
```

## Download FastQC in the environment
source: https://anaconda.org/bioconda/fastqc
```
conda install -c bioconda fastqc
conda install -c bioconda/label/broken fastqc
```
## Run FastQC on all fastq files and specify output folder
Once all the fastq files are downloaded, run the following command when inside your readqc folder.
```
cd readqc
mkdir fastqc_outputs
fastqc *.fastq --outdir=./fastqc_outputs
```

## Deactivate Environment
```
source deactivate
```

## Terminate Instance
On EC2 terminal tab, terminate the instance. Check online through aws to make sure the instance has been fully terminated.
```
aegea terminate kalani-cc4
```
