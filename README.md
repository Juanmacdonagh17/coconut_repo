***Coconut is currently being developed.
Be aware of bugs, especially tyranids.***

---

![alt text](https://github.com/Juanmacdonagh17/coconut_repo/blob/main/logo/logo2.svg)

# Coconut: 
## A suite for transcripts, codon usage and protein structure analysis.
---

### Introduction 

Coconut is an easy-to-run codon analysis suite. It's fully coded on C, providing functions for analyzing FASTA or multi FASTA files or fetching transcript sequences from Ensemlb, and calculating parametires like Codon Usage, Relative Synonimus Codon Usage and %MinMax.  

---
### Installation 

The software is easily complied via [GCC](https://gcc.gnu.org/). To install, you need to have the curl.h library.
We provide a shell script for quick installation:

The script checks for GCC, wget and make and curl. If they are installed, it simply compiles coconut. If not, it offers to install miniconda and compiles coconut.
```bash
git clone https://github.com/Juanmacdonagh17/coconut_repo
cd coconut_repo
chmod +x coco_install.sh
./coco_install.sh
```

If you prefer to install it by yourself, you can install [Miniconda3](https://docs.anaconda.com/miniconda/) (or simply just install [curl](https://github.com/curl/curl), although you might need some extra dependencies).
After that, just run the following command:

```bash
git clone https://github.com/Juanmacdonagh17/coconut_repo
cd coconut_repo
gcc -o coconut coconut.c -I/home/path/to/miniconda3/include -L/home/path/to/miniconda3/lib -lcurl

```
Replace /path/to with the path  where miniconda was installed.

Coconut also counts with a Docker implementation:
Be sure to have [Docker](https://www.docker.com/) installed on your computer, then follow this instructions:

```bash
git clone https://github.com/Juanmacdonagh17/coconut_repo
cd coconut_repo
docker build -t coconut-app .
```
After that, you can run Coconut via docker like this:
```bash
docker run coconut-app -help
```
Mounting a volume:
```bash
docker run -v /home/user/coconut_output:/home/cocouser/output your_image_name -help

```
---

### Running Coconut 

To run Coconut, after installing, just execute ./coconut in the command line.
Using the -help flag will provide a list of all available commands, but here is a list of commands that Coconut can use (so far):

+ -help: shows all the arguments and options available.
+ -slinet: optional flag, used for not showing the outputs while running Coconut.
+ -cu: calculates Codon Usage for an input fasta or multifasta or a fetched ENST.
+ -rscu: calculates Relative Synonimus Codon Usage for an input fasta or multifasta or a fetched ENST.
+ -all: calculates both CU and RSCU for an input fasta or multifasta or a fetched ENST.

+ -fetch: fetches a transcript (specifically the main CDS) from Ensembl, using the transcript identifier.
+ -fetchfile: same as fetch, but using a list of identifiers from a file, one per line.
+ -multi: fetches all versions of a transcript, using ENSG instead of ENST.
+ -gene: fetches all versions of a transcript using a species and gene symbol as inputs. 
+ -uniprot: fetches the main version of a transcript using the UniProt ID as an entry point.
+ -af: fetches the pLDDT from the AF model, and it creates a naive structural assignation based on pLDDT and contacts from the PDB AF file.
+ -slice_domains: used to easily compare CU and RSCU between domains. It uses a file with "slice instructions" (find examples below, works with local files and with fetch).

+ -minmax: calculates the %MinMax for an input fasta or multifasta. If not indicated otherwise, the window of comparison is 18 codons, but this can be changed. For this, a table of optimal use codons is needed. The table for _Homo Sapiens_ is provided.

---
### Examples

There are some examples files provided in this repository: 

To run Coconut with an input and output file you can run:

```bash
./coconut -cu multi_cds.fasta multi_cds.csv
```
If you want to fetch sequences using their Ensembl transcript ID from a file:

```bash
./coconut -all -fetchfile ensts.txt ensts.csv
```
And for a single file: 

```bash
./coconut -rscu -silent -fetch ENST00000361390 ENST00000361390.csv
```
Fetching an UniProt and the AlphaFold model (this also creates a file with pLDDT, contacts and a naive domain classification):

```bash
./coconut -all -af -uniprot P05067 aa.csv
```

To perform cu/rscu over the AF domain regions:

```bash
./coconut -cu -af_regions -uniprot O95905 O95905.csv
```

To use the slice function (an example slice.csv is provided to show how the file needs to look like):

```bash
./coconut -cu -slice_domains  example.fasta slice.csv example.csv
```
Fetching all versions of a transcript using a gene symbol: 

```bash
./coconut -cu -gene homo_sapiens INS insutlin_output.csv
```
To perform %MinMax calculations (no output name is required here):

```bash
./coconut -minmax Usage-num3.man multi_cds.fasta
```
And if you want to use a different window of codons:

```bash
./coconut -minmax  Usage-num3.man 10 multi_cds.fasta
```

If you want to perform the RRT over the sequence:

```bash
./coconut -rtt -minmax Usage-num3.man ENST00000372979
```
---

### Progress: 

| Feature                                     | Available |
|---------------------------------------------|-----------|
| Read FASTA                                  | ✅        |
| Read multi-FASTA                            | ✅        |
| Calculate CU                                | ✅        |
| Calculate RSCU                              | ✅        |
| Fetch ENST from Ensembl                     | ✅        |
| Fetch multiple ENSTs from list (Ensembl)    | ✅        |
| Slice domains from FASTA/multi-FASTA        | ✅        |
| Calculate MinMax                            | ✅        |
| Docker image                                | ✅        |
| Fetch multiple versions of a transcript     | ✅        |
| Fetch using gene and species symbol         | ✅        |
| Slice domains from fetch                    | ✅        |
| Slice multiple versions from fetch          | ✅        |
| UniProt as entry point for request          | ✅        |
| AFDB integration for structural data        | ✅        |
| AF structure classification                 | ✅        |
| Calculate Codon Adaptation Index            | ✅        |
| Perform RRT over the %MinMax                | ✅        |
| Perform CU/RSCU over AF domain classifiation| ✅        |
| Fetch domains from InterPro or CATH         | ❔        |
| ~~NW alignment method~~                     | ❌        |
| ~~MobiDB integration for disorder regions~~ | ❌        |
| ~~Incorrect length error~~                  | ❌        |
| ~~Optional flags (no ATG, no STOP, etc)~~   | ❌        |
                                                      
---



