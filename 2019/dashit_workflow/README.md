# DASH Introduction
Got junk sequences in your sequencing fastqs? Many RNA-seq experiments are riddled with repetitive sequences which are of no value to the scientific question you are asking. Some examples include ribosomal RNA and hemoglobin. DASH is a CRISPR-cas9, NGS technology which depletes abundant sequences, increasing coverage of your sequences of interest.

In the lab, you take final NGS libraries ready for sequencing, and combine them with cas9-gRNAs which cut up your repetitive sequences into small fragments. You then size select and amplify your remaining DNA, which should contain most or only your sequences of interest. Your library can then be sequenced. For more details, see our protocol on protocols.io: dx.doi.org/10.17504/protocols.io.y9gfz3w

# What is DASHit?
In order for DASH to be effective, you want to pick the fewest guides which hit the most of your repetitive sequences. The best way to know the most repetitive sequences in your sequencing is by designing guides which hit your actual reads. David has implemented a greedy algorithm to solve this set cover problem, and created a software package called DASHit to perform automated guide design for DASH experiments.

DASHit uses your reads from a preliminary, low-depth sequencing of your samples, and identifies crispr-cas9 cut sites. It has a filtering function for the guides it identifies, which allows you to specify on and off target sequences, and filter based on GC content and secondary structure. It then optimizes your guide list by giving you a set number of guides which hit the most reads.

DASHit also has functionality to score guides against a fasta. You can use this to test how many of your reads would be hit by a given guide set.

# Cupcakes 05/16/19 - let's get started!
Today, we will be using a set of preliminary reads on mouse tissue to design a set of 96 DASH guides. I have written a wrapper script which encompasses most of this work, which will be made available on the DASHit github repository, and remain up-to-date as a reference for anyone wishing to run this pipeline.

## Launching an AWS instance
The czbiohub-specops image has dashit installed on it. However, that is an outdated version, so we will be running dashit on a clean instance, with an image I made containing seqtk, bedtools and cutadapt.

Replace the instance name with your name.

```
 aegea launch --iam-role S3fromEC2 --ami-tags Name=dashit-image-v3 -t m4.large --no-dns dashit-yourname --duration-hours 2
```

### Login to your AWS instance
```
aegea ssh ubuntu@dashit-yourname
```

Change directory into volume /mnt/data
```
cd /mnt/data/
```

## Install DASHit (SKIP SINCE USING IMAGE)
Installing DASHit requires python3 and GO.

```
aws s3 cp s3://alyden-bucket/Mouse_DASH/cupcakes/DASHit_installation.sh . --region "us-east-2"
```

You will be prompted to read the license agreement. Click q when you are done, and type I agree, and hit enter.

```
bash DASHit_installation.sh
```

Here are the contents of the script:
```
wget https://dl.google.com/go/go1.12.5.linux-amd64.tar.gz
sudo tar -C /usr/local -xzf go1.12.5.linux-amd64.tar.gz
export PATH=$PATH:/usr/local/go/bin
python3 -m venv ~/.virtualenvs/dashit
source ~/.virtualenvs/dashit/bin/activate
git clone https://github.com/czbiohub/guide_design_tools.git
cd guide_design_tools
make install
ln -s /mnt/data/guide_design_tools/vendor/special_ops_crispr_tools/offtarget/offtarget /home/ubuntu/bin/offtarget
ln -s /mnt/data/guide_design_tools/vendor/special_ops_crispr_tools/crispr_sites/crispr_sites /home/ubuntu/bin/crispr_sites
```

## Enter DASHit virtual environment
```
source ~/.virtualenvs/dashit/bin/activate
```

## Copy files from AWS

Make a new directory and change into it
```
mkdir DASH
cd DASH
```

Sync files needed for this from AWS. This will take a few minutes.
```
aws s3 sync s3://alyden-bucket/Mouse_DASH/cupcakes/ . --region "us-east-2"
```

`ls` to take a look at what you downloaded! There should be a series of paired-end reads, two genome fastas and a bed file.

## Preparing your reads to run DASHit
## SUBSAMPLE

**Subsample files**: if your files are large, we recommend randomly subsampling your files to 100k reads. This can be done on a fastq or fastq.gz using `seqtk`. Use the same seed to maintain paired information.

### Install seqtk (SKIP SINCE USING IMAGE)

```
sudo dpkg --configure -a
sudo apt install seqtk
```

**Example commands**
```
seqtk sample -s100 input_R1.fastq 100000 > sub100k_input_R1.fastq
```
```
seqtk sample -s100 input_R2.fastq 100000 > sub100k_input_R2.fastq
```
### Subsample all of your reads to 10000 in a for loop
```
for i in GI01_*fastq.gz; do seqtk sample -s100 $i 10000 > sub10k_${i:0:-3}; done
```

## TRIM ADAPTORS
This step is crucial for maintaining library integrity! The last thing we want is to design a guide which cuts your adaptor sequence in it, which would render your library unsequencable. To avoid this, trim all adaptors off your reads using your favorite bioinformatics tool. We do this using `cutadapt`, but other programs such as `trimmomatic` or `fastp` also work well.

### Install cutadapt (SKIP SINCE USING IMAGE):

```
pip install cutadapt
```
or
```
sudo apt install python-cutadapt
```

**Example command**
```
cutadapt --report=minimal -j 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o cut-sub10k_input_R1.fastq -p cut-sub10k_input_R2.fastq sub10k_input_R1.fastq sub100k_input_R2.fastq
```

### Cut adaptors off all of your reads in a for loop
```
for i in sub10k*R1*; do cutadapt -j 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o cut-$i -p cut-${i:0:-11}2_001.fastq $i ${i:0:-11}2_001.fastq -m 40; done
```

## CONVERT TO FASTA
**Convert from fastq to fasta if needed**: you can convert from fastq to fasta in many ways; we use `seqtk` for this as well.

**Example commands**
```
seqtk seq -A  cut-sub100k_input_R1.fastq > cut-sub100k_input_R1.fasta
```
```
seqtk seq -A  cut-sub100k_input_R2.fastq > cut-sub100k_input_R2.fasta
```
### Convert all reads from fastq to fasta in a for loop
```
for i in cut-sub10k_GI01_*fastq; do seqtk seq -A $i > ${i:0:-1}a; done
```

Let's take a look and confirm they are fastas with trimmed adaptors:
```
head *fasta
```

## Running crispr_sites to identify all PAM sites in your reads
```
crispr_sites -h
```
crispr_sites will take your reads from stdin and output a list of guides. If you use -r flag (for reads), it will track which read it hit.

```
cat cut-sub10k_GI01_*fasta | crispr_sites -r > mouse_reads_crispr_sites.txt
```

Take a look at the crispr_sites file!

## Filtering guides with dashit-filter (optional)
If you want to use these guides as is, you can run straight into optimize guides, to identify the guides which hit the most reads, using the above crispr_sites file. However, most people want to filter their guides for on or off target activity.

Let's filter by both on and off target activity.

**On-target**: I want to hit the ribosomal regions of the mouse.

**Off-target**: I do not want to hit the rest of the mouse genome, because I am interested in mouse transcriptome, and I also do not want to hit the *Candida albicans* genome, because I am trying to quantify infection.

### Preparing your on and off target files

You can use `bedtools` to help generate your on-target and off-target fastas.

## Install bedtools (SKIP SINCE USING IMAGE)
```
sudo apt install bedtools
```
Note: having trouble with sudo sometimes on instances because username not in /etc/host, works well otherwise...try some of these options if it is not working https://bedtools.readthedocs.io/en/latest/content/installation.html

**On-target files**

If you have a bed file annotating on-target regions, you can use `bedtools getfasta` to generate an on-target fasta. I created a bed file from UCSC Genome Browser, pulling all regions which were annotated as ribosomal in the mm10.fa genome. I then downloaded the mm10.fa.

Using these two files, bedtools will pull out the on-target regions of the mm10 genome.
```
bedtools getfasta -fi mm10.fa -bed mm10_rRNA.bed -fo on-target_mm10_rRNA.fa
```

**Off-target files**
You can use the same bed file to mask on-target regions of a fasta to generate an off-target fasta, using `bedtools maskfasta`, with the `-mc -` option. This will replace all regions in the bed file with a `-` instead of a nucleotide. Since I want to preserve other elements of the mouse genome, I will mask everything that is not rRNA, as annotated by my bed file. [This step will take 1-3 minutes]

```
bedtools maskfasta -bed mm10_rRNA.bed -fi mm10.fa -fo off-target_mm10.fa -mc -
```

Quick check:
```
grep "-" off-target_mm10.fa | head
```

Since I am also hoping to avoid hitting the *C. albicans* genome with my guides, I will add that to my off-target file. [This step will take 1-3 minutes]
```
cat off-target_mm10.fa CalbicansSC5314_genome_allNs.fasta > off-target_mm10_SC5314.fa
```

*Note*: Both fastas should contain no letters other than A, G, C, T or N.

###Running crispr_sites on your on and off target fastas
We will now run crispr_sites (like before), without the -r flag since these are not reads. This may take 4-5 minutes.

```
cat on-target_mm10_rRNA.fa | crispr_sites > on_target_crispr_sites.txt
cat off-target_mm10_SC5314.fa | crispr_sites > off_target_crispr_sites.txt
```

###Running dashit-reads-filter to include your on-target guides and exclude off-target guides

We can now use our on-target crispr sites file, our off-target crispr sites file and our reads crispr sites file to run dashit-reads-filter. Run the following command to see all the options for filtering, including structural and GC content flags.

```
dashit-reads-filter -h
```

usage: dashit-reads-filter [-h] [--filtered_explanation FILTERED_EXPLANATION]
                           [--offtarget OFFTARGET]
                           [--offtarget_radius OFFTARGET_RADIUS]
                           [--ontarget ONTARGET]
                           [--ontarget_radius ONTARGET_RADIUS]
                           [--gc_freq_min GC_FREQ_MIN]
                           [--gc_freq_max GC_FREQ_MAX]
                           [--homopolymer HOMOPOLYMER]
                           [--dinucleotide_repeats DINUCLEOTIDE_REPEATS]
                           [--hairpin_min_inner HAIRPIN_MIN_INNER]
                           [--hairpin_min_outer HAIRPIN_MIN_OUTER]
                           input

Filter guides in a sites-to-reads file based on offtargets and quality

positional arguments:
  input                 input sites-to-reads file to filter. Generated by
                        crispr_sites -r

optional arguments:
  -h, --help            show this help message and exit
  --filtered_explanation FILTERED_EXPLANATION
                        output file listing which guides were disqualified and
                        why. CSV format.

offtarget filtering:
  options to filter offtargets

  --offtarget OFFTARGET
                        File containing off target CRISPR sites, as generated
                        by crispr_sites
  --offtarget_radius OFFTARGET_RADIUS
                        Radius used for matching an off target. Specify this
                        as L_M_N which means remove a guide for hitting an off
                        target if L, M, N nucleotides in the first 5, 10 and
                        20 positions of the guide, respectively, match the off
                        target. e.g., 5_10_20 to require perfect matches;
                        5_9_18 to allow up to one mismatch in positions 6-10
                        positions and to allow up to 2 mismatches in the last
                        10 positions

ontarget filtering:
  options to filter ontargets

  --ontarget ONTARGET   File containing ontarget CRISPR sites, as generated by
                        crispr_sites
  --ontarget_radius ONTARGET_RADIUS
                        Radius used for matching ontargets. Same format as
                        --offtarget_radius.

quality filtering:
  options for how guides are filtered for poor structure reasons

  --gc_freq_min GC_FREQ_MIN
                        filter guide if # of Gs or Cs is strictly less than
                        this number
  --gc_freq_max GC_FREQ_MAX
                        filter guide if # of Gs or Cs is strictly greater than
                        this number
  --homopolymer HOMOPOLYMER
                        filter guide if strictly more than this number of a
                        single consecutive nucleotide appears, e.g., AAAAA
  --dinucleotide_repeats DINUCLEOTIDE_REPEATS
                        filter guide if strictly more than this number of a
                        single dinucleotide repeats occur, e.g. ATATAT
  --hairpin_min_inner HAIRPIN_MIN_INNER
                        filter guide if a hairpin occurs with >=this inner
                        hairpin spacing, e.g., oooooIIIooooo, where the o are
                        reverse complements and III is the inner hairpin
                        spacing
  --hairpin_min_outer HAIRPIN_MIN_OUTER
                        filter guide if a hairpin occurs with >=this outer
                        hairpin width, e.g., oooooIIIooooo, where the o are
                        reverse complements and ooooo is the outer hairpin


On default settings, a site will be removed if any of the following are true:

1. G/C frequency too high (> 15/20) or too low (< 5/20)
2. Homopolymer: more than 5 consecutive repeated nucleotides
3. Dinucleotide repeats: the same two nucelotides alternate for > 3
   repeats
4. Hairpin: complementary subsequences near the start and end of a
   site can bind, causing a hairpin

We will require perfect matches to define on and off target activity.

```
dashit-reads-filter mouse_reads_crispr_sites.txt --offtarget off_target_crispr_sites.txt --offtarget_radius 5_10_20 --ontarget on_target_crispr_sites.txt --ontarget_radius 5_10_20 --filtered_explanation why_final_mouse_rRNA_guides.csv > mouse_rRNA_final_crispr_sites.txt
```

## Optimize guides
Now that we have a filtered guide list, we can use this to pick the fewest guides which will hit the greatest number of reads. We run optimize guides, which allows us to specify the number of desired guides, as well as how many times a read has to be cut before it is considered 'hit'. For us, we will design 500 and say each read must only be hit 1 time.

```
optimize_guides mouse_rRNA_final_crispr_sites.txt 500 1 > mouse_rRNA_DASH_500_guides.csv
```

Take a look at this file! How much can we DASH with 500 guides? 250? 100?

The file will contain five components.
1. Guide (number)
2. Site (the sequence)
3. Number of reads covered by site
4. cumulative number of reads covered
5. cumulative percent of reads covered 

```
less mouse_rRNA_DASH_500_guides.csv
```

Let's take the first 96 guides and score them against our original files.
```
head -98 mouse_rRNA_DASH_500_guides.csv > mouse_rRNA_DASH_96_guides.csv
```

## Score guides
We can check *in silico* see how each file would be DASHed if we were to use this guide set.

**Example command**
```
score_guides mouse_rRNA_DASH_96_guides.csv sub10k_GI01_7_S3_R2_001.fasta
```

Run score_guides in a for loop to check out all of the files:

```
for i in cut-sub10k_GI01_*fasta; do score_guides mouse_rRNA_DASH_96_guides.csv $i >> mouse_rRNA_final96_guides_scored.txt; done;
```

Take a look at the `mouse_rRNA_final96_guides_scored.txt` file!

If you want, you can format this txt file into a CSV using a custom Python script:
```
python3 ../guide_design_tools/dashit/contrib/score_guides_scripts/DASH_csv_format.py mouse_rRNA_final96_guides_scored.txt
```
Your formatted CSV will contain five columns
1. Guide library
2. Your filename
3. Total Reads DASHed (hit by score_guides)
4. Total Reads in sample
5. Percent DASHed

```
less mouse_rRNA_final96_guides_scored.csv
```

### Running the design_guides_wrapper script

This script will run crispr_sites, optimize_guides, score_guides, and the formatting Python script. Once you have created your on-target and off-target fastas, you can input the following into the wrapper:

Use bash to run the script and follow with 3 arguments
1. the file prefix for all of your read fastas
2. the file name of your on-target fasta
3. the file name of your off-target fasta

```
bash design_guides_wrapper.sh cut-filt-85-98-90_sub100k_GI01 ontarget_mm10-rRNA_region_UPPERcase.fa masked_mm10_SC5314_genomes_UPPERcase.fa
```

### Your final guides/optimize_guides output file
Your final guides CSV (`final_guides.csv`) will contain five components.
1. Guide (number)
2. Site (the sequence)
3. Number of reads covered by site
4. cumulative number of reads covered
5. cumulative percent of reads covered 

You should use this file to generate an elbow plot, with guide sequence on the x-axis and number of reads covered by site/total reads hit on the y-axis.

### Your formatted score_guides output file

We will then run optimize_guides to determine the number of reads hit by each guide, and run score_guides on your original reads with the full guide set to see how much is DASHable. The optimize_guides output can be used to make an elbow plot to determine the optimal number.

Your formatted CSV (`final_guides_scored.csv`) will contain five columns
1. Guide library
2. Your filename
3. Total Reads DASHed (hit by score_guides)
4. Total Reads in sample
5. Percent DASHed

This file will tell you how much of your library would be DASHed by the full guide set. Use the elbow plot to determine how many guides would be most effective and rerun score_guides.

### Other files created
Files with `crispr_sites.txt` will contain guide sequences in your reads, on-target or off-target fastas which will be compatible with optimize_guides. A txt file called `final_guides_scored.txt`  is generated containing score_guides output before it is formatted into a CSV.
