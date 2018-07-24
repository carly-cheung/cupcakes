# Launch EC2 Instance with Bowtie2
Bowtie2 is installed on the Special Ops image in the conda environment _py2_.
`aegea launch --iam-role S3fromEC2 --ami-tags Name=czbiohub-specops -t m5.4xlarge shao-cc`

# Connect to Instance and Download Data from S3
```bash
aegea ssh ubuntu@shao-cc
cd /mnt/data
aws s3 cp s3://czbiohub-seqbot/fastqs/180615_NB501961_0126_AHWV5YAFXX/rawdata/RR011_006_iso_rna1_H52678_S6_R1_001.fastq.gz
aws s3 cp s3://czbiohub-seqbot/fastqs/180615_NB501961_0126_AHWV5YAFXX/rawdata/RR011_006_iso_rna1_H52678_S6_R2_001.fastq.gz
```
# Upload Reference Genome from Local
- reference genome from NCBI contained in this repo as "Staph\_aureus\_reference.fasta"
In your local terminal, navigate to the directory containing the file.
Upload the file to your instance, e.g. `aegea scp ubuntu@shao-cc:/mnt/data`
Go back to your instance.

# Tmux or Screen
We want to start Bowtie2 in a tmux/screen session since this lets us multiplex the screen.

Quick-launch tmux session: `tmux`
Named tmux session: `tmux new -s bowtie`

_Note:_ This tutorial will be done using tmux key bindings. Google "gnu screen cheat sheet" to find the corresponding screen key bindings.

# Activate Conda Environment
Bowtie2 is installed in the conda environment _py2_ (SRST2, which depends on Bowtie2 and Python 2, is installed here).

`source activate py2`

# Using Bowtie2
## Get Manual Without Going to Browser
For most bioinformatics tools, either typing the command or adding the flag `--help` or `-h` will bring up a usage manual.

`bowtie2 --help` or `bowtie2 -h` or `bowtie2`

## Taking Advantage of Screen/Tmux
Once we have the manual, we'd like to be able to refer to it as we're typing our bowtie2 command. To do this, we'll split our screen vertically.

```bash
^B      # All tmux commands are preceded by ctrl + b
%       # key binding to vertically split screen
```
In order to move around on the window, we have to enter copy mode.
```bash
^B
[
# Use arrow keys to move up and down. Hold Fn to move by page.
```
Now switch to the new pane, where we will actually run Bowtie2.

```bash
^B
->      # Use arrow keys to move between panes
```
Google "tmux cheat sheet" to find the key bindings for more functions.

## Build Index
Bowtie2 requires an index built from our reference.
- To see usage: `bowtie2-build -h`
- To build our index, type `bowtie2-build Staph_aureus_reference.fasta S_aur_ref`.

## Aligning Reads to Reference
You should refer to the manual to see what options you need. Below is what a default Bowtie2 command looks like.

```bash
bowtie2 -x S_aur_ref -1 RR011_006_iso_rna1_H52678_S6_R1_001.fastq.gz -2 RR011_006_iso_rna1_H52678_S6_R2_001.fastq.gz -S sample1_S_aur.sam
```

- Since we're working in an EC2 instance, we want to take advantage of the multithreading options: `-p 16`
- We'll also only look at the alignment information for reads that aligned, and not bother with the unaligned reads: `--no-unal`
- For this tutorial, we'll use a fast end-to-end alignment. We might want more local or more sensitive alignment depending on our goals: `--very-fast`
- Bowtie2 outputs overall alignment statistics to stderr, so let's redirect that to a file: `2> sample1_S_aur.out`

```bash
bowtie2 -x S_aur_ref -1 RR011_006_iso_rna1_H52678_S6_R1_001.fastq.gz -2 RR011_006_iso_rna1_H52678_S6_R2_001.fastq.gz -S sample1_S_aur.sam -p 16 --no-unal --very-fast 2> sample1_S_aur.out
``` 

## Observing Bowtie2
Now that we're done with the manual, we can use that pane to observe our process instead.
In the other pane:

```bash
q       #exits copy mode
sudo apt install htop
htop        # Can also just use top, which is already installed
```

## Reading Overall Alignment Statistics
```bash
cat sample1_S_aur.out
```
or
```bash
less sample1_S_aur.out
q       # to quit
```

# Primer on SAM format
Good quick guide to SAM/BAM [here](https://training.h3abionet.org/postgraduate_workshop_2014/wp-content/uploads/2014/04/H3ABioNet_2014_NGS_8_SamFormat.pdf).
- Can't load SAM/BAM to view in Geneious without reference FASTA

# Terminate Instance
`aegea terminate shao-cc`
