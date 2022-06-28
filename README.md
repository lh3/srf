## Getting Started

```sh
# compile srf
git clone https://github.com/lh3/srf
cd srf && make

# count high-occurrence k-mers with KMC
ls *.fastq.gz > fofn.txt && mkdir -p tmp_dir
kmc -fq -k151 -t16 -ci100 -cs100000 @fofn.txt count.kmc tmp_dir
kmc_dump count.kmc count.txt

# assemble satellite DNA
./srf -p prefix count.txt > srf.fa
```

## Introduction

Satellite Repeat Finder, or srf in brief, assembles motifs in [satellite
DNA][satdna] that are tandemly repeated many times in the genome.

[satdna]: https://en.wikipedia.org/wiki/Satellite_DNA
