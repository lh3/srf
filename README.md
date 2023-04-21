## Getting Started

```sh
# compile srf
git clone https://github.com/lh3/srf
cd srf && make

# count high-occurrence k-mers in HiFi reads with KMC
ls *.fastq.gz > fofn.txt && mkdir -p tmp_dir
# set -ci to approximately 10-fold over the average coverage
kmc -fq -k151 -t16 -ci100 -cs100000 @fofn.txt count.kmc tmp_dir
kmc_dump count.kmc count.txt

# assemble satellite DNA
./srf -p prefix count.txt > srf.fa

# count high-occurrence k-mers in a HiFi or duplex assembly
kmc -fm -k151 -t16 -ci20 -cs100000 ref.fa count.kmc tmp_dir
kmc_dump count.kmc count.txt

# analyze
minimap2 -c -N1000000 -f1000 -r100,100 <(./srfutils.js enlong srf.fa) ctg.fa > srf-aln.paf
./srfutils.js paf2bed srf-aln.paf > srf-aln.bed   # filter and extract non-overlapping regions
./srfutils.js bed2abun srf-aln.bed > srf-aln.len  # calculate abundance of each srf contig
```

## Introduction

Satellite Repeat Finder, or SRF in brief, assembles motifs in [satellite
DNA][satdna] that are tandemly repeated many times in the genome. It takes
short reads, accurate long reads or high-quality contigs as input and reports
the consensus of each repeat unit. SRF can identify satellite repeats that are
often missed in de novo assembly. For species enriched with high-order repeats
(HORs), it tends to find HORs instead of the minimal repeat unit. SRF may also
find truly circular genomes such as mitochondial or chloroplastic genomes if
their abundance is high.

## Usage

### Input

SRF works best with phased telomere-to-telomere assemblies and may work with
trio hifiasm assemblies. If you worry about non-uniform assembly quality, you
may use accurate long reads (e.g. PacBio HiFi or nanopore duplex)
as input. The down side is that you would not be able to obtain phased results.
SRF also works with short reads at reduced power.

### Generating satellite contigs

To use SRF, first count high-occurrence k-mers. We use [KMC][kmc] as an example:
```sh
kmc -fq -k151 -t16 -ci100 -cs100000 HiFi-reads.fq.gz count.kmc tmp_dir
kmc_dump count.kmc count.txt
```
We recommend to use a large `-k` for long k-mers and set a high `-ci` to skip
most k-mers in the unique regions of the genome. The proper choice of these
parameters varies with the input data. For short reads, `-k` should be much
shorter than the read length to reach enough k-mer coverage. For HiFi reads, we
set `-ci` to 10 times the average read coverage in [our preprint][pre]. For
contigs, we used `-ci20`. Highly abundant satellite repeats are usually
insensitive to the KMC parameters as long as k>100 but less abundant repeats
are not as stable.

After k-mer counting, run SRF to get contigs:
```sh
./srf -p prefix count.txt > srf.fa
```
This step usually takes a couple of minutes. The output looks like:
```txt
>prefix#circ1-2126 min=34936,max=100000,avg=94336
AGACCGTGGGGATGCTGGCGAAGCT...
>prefix#circ2-2214 min=6668,max=76009,avg=20245
GGCTTCTTCCCTTGAGCTCTGCAGC...
...
```
The first FASTA record tells us the first contig is 2126bp in length. The
minimum, maximum and average k-mer occurrences on the contig are 34936, 100000
and 94336, respectively.

### Analysis

To estimate the abundance of each satellite repeat, we may map SRF contigs
back to the source sequences (requiring the k8 javascript engine):
```sh
./srfutils.js enlong srf.fa > srf.enlong.fa  # enlong short contigs for mapping
minimap2 -c -N1000000 -f1000 -r100,100 -t16 srf.enlong.fa HiFi-reads.fa > srf.paf
./srfutils.js paf2bed srf.paf > srf.bed  # generate non-redundant mapping regions
```
In `srf.paf`, two SRF contigs may be mapped to the same locus. The `paf2bed`
command of `srfutils.js` resolves such redundancy with a simple heuristic. We
can estimate the abundance with
```sh
./srfutils.js bed2abun srf.bed
```

SRF may report HORs composed of similar monomers. These may not be stable HORs.
You may use [TRF](https://tandem.bu.edu/trf/trf.html) or
[TRF-mod](https://github.com/lh3/TRF-mod) to decompose such HORs to monomers.
SRF may also report multiple similar sequences for the same type of a repeat.
You may run
```sh
minimap2 -c -N1000 <(./srfutils.js enlong -d srf.fa) srf.fa
```
and filter SRF sequences with your own criteria.

### Citation

If you use SRF, please consider to cite the following preprint

> Zhang Y, Chu J, Cheng H and Li H (2023) De novo reconstruction of satellite
> repeat units from sequence data. [arXiv:2304.09729][pre]

[satdna]: https://en.wikipedia.org/wiki/Satellite_DNA
[kmc]: https://github.com/refresh-bio/KMC
[pre]: https://arxiv.org/abs/2304.09729
