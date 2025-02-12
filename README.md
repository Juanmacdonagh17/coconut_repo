***Coconut is currently being developed.
Be aware of bugs***

---

![alt text](https://github.com/Juanmacdonagh17/coconut_repo/blob/main/logo/logo2.svg)

# Coconut: 
## A suite for transcripts, codon usage and protein structure analysis.
---

### Introduction 

Coconut is an easy-to-run codon bias and protein structure analysis suite. It's fully coded on C, providing functions for analyzing FASTA or multi FASTA files or fetching transcript sequences from Ensemlb or UniProt, and calculating parameters like Codon Usage, Relative Synonymous Codon Usage, %MinMax and simple protein structure analysis.  

---
### Biological context

According to the central dogma of molecular biology, every protein comes from an RNA sequence, that, in turn, comes from a DNA sequence. 

It is already established that the function of proteins is strictly correlated with the - or lack of - structure and that protein structure can be determined by the contact of the atoms that form the aminoacidic chain. 
Recently, it has been proven that not only the aminoacidic determines the structure - or ensemble of structures - but that the use of codons at the messenger RNA level can also affect this, as they regulate the speed  of translation within the ribosome, and can give the protein enough time to fold into separate domains (for more context, see references (1, 2)).

Even though many online tools and packages for different programming languages can calculate different parameters to study these correlations, there's still one single suite of tools that can perform this different analysis in a fast and easy to run way. With this C pipeline, we intend to provide the bioinformatics and molecular biology community with a single tool that can perform reliable analysis for these cases.

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
+ -rscu: calculates Relative Synonimus Codon Usage for an input fasta or multifasta or a fetched ENST. (For more context on RSCU and CAI, check references (3).
+ -all: calculates both CU and RSCU for an input fasta or multifasta or a fetched ENST.

+-cai: calculates the Codon Adaptation Index, using a codon reference table with frequencies. (For more context on RSCU and CAI, check references (3).

+ -fetch: fetches a transcript (specifically the main CDS) from Ensembl, using the transcript identifier.
+ -fetchfile: same as fetch, but using a list of identifiers from a file, one per line.
  
+ -multi: fetches all versions of a transcript, using ENSG instead of ENST.
+ -gene: fetches all versions of a transcript using a species and gene symbol as inputs.
  
+ -uniprot: fetches the main version of a transcript using the UniProt ID as an entry point.
+ -af: fetches the pLDDT from the AF model, and it creates a naive structural assignation based on pLDDT and contacts from the PDB AF file. (For more context on how the regions are assigned, check references (5, 6).
+ af_regions: fetches AlphaFold pLDDT data and slice regions to a file using UniProt ID

+ -slice_domains: used to easily compare CU and RSCU between domains. It uses a file with "slice instructions" (find examples below, works with local files and with fetch).
  

+ -minmax: calculates the %MinMax for an input fasta or multifasta. If not indicated otherwise, the window of comparison is 18 codons, but this can be changed. For this, a table of optimal use codons is needed. The table for _Homo Sapiens_ is provided. (For more context on how the %MinMax is calculated, check references (4).
+ +- rrt: calculate the Random Reverse Translations for the %MinMax input sequence (1000 iterations). (For more context on how the %MinMax is calculated, check references (4).


Aditional codon usage tables can be found online, for example: https://www.kazusa.or.jp/codon/ 
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
| Calculate CAI                               | ✅        |
| Fetch domains from InterPro or CATH         | ❔        |
                                                      
---
### References: 

1. The Effects of Codon Usage on Protein Structure and Folding
McKenze J. Moss, Laura M. Chamness, and Patricia L. Clark
2. Synonymous but Not Silent: The Codon Usage Code for Gene Expression and Protein Folding
Yi Liu, Qian Yang, and Fangzhou Zhao
3. The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications
Paul M. Sharp, Wen-Hsiung Li
4. %MinMax: A versatile tool for calculating and comparing synonymous codon usage and its impact on protein folding
Anabel Rodriguez, Gabriel Wright, Scott Emrich, Patricia L. Clark
5. Accurate structure prediction of biomolecular interactions with AlphaFold 3
Josh Abramson et. al.
6. Afflecto: A Web Server to Generate Conformational Ensembles of Flexible Proteins from Alphafold Models
Matyas Pajkos, Ilinka Clerc, Christophe Zanon, Pau Bernado, Juan Cortes

---


