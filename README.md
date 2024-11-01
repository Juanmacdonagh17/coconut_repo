***Coconut is currently being developed.
Be aware of bugs, especially tyranids.***

---

![alt text](https://github.com/Juanmacdonagh17/coconut_repo/blob/main/logo/logo2.svg)

# Coconut: 
## A suite for transcripts & codon usage calculations
---

### Introduction 

Coconut is an easy-to-run codon analysis suite. It's fully coded on C, providing functions for analyzing FASTA or multi FASTA files or fetching transcript sequences from Ensemlb, and calculating parametires like Codon Usage, Relative Synonimus Codon Usage and %MinMax.  

---
### Instalation 

The software is easily installed via [GCC](https://gcc.gnu.org/). To install, you need to have the curl.h library.
The simplest, if you don't have it already, is installing [Miniconda3](https://docs.anaconda.com/miniconda/).

After that, just run the following command:


```bash
git clone https://github.com/Juanmacdonagh17/coconut_repo
cd coconut_repo
gcc -o coconut coconut.c -I/home/path/to/miniconda3/include -L/home/path/to/miniconda3/lib -lcurl

```
Replace /path/to with the folder where Conda was installed.

Coconut also counts with a Docker implementation:

Be sure to have [Docker](https://www.docker.com/) installed on your computer, then follow this instructions:

```bash
git clone https://github.com/Juanmacdonagh17/coconut_repo
cd coconut_repo
docker build -t coconut-app .
```
After that, you can run Coconut via docker like this:
```bash
docker run -v /usr/src/coconut_repo ARGUMENTS
```
---

### Running Coconut 

To run Coconut, after installing, just execute ./coconut in the command line.
Using the -help flag will provide a list of all available commands, but here is a list of commands that Coconut can use (so far):

+ -help: shows all the arguments and options available
+ -slinet: optional flag, used for not showing the outputs while running Coconut.
+ -cu: calculates Codon Usage for an input fasta or multifasta or a fetched ENST.
+ -rscu: calculates Relative Synonimus Codon Usage for an input fasta or multifasta or a fetched ENST.
+ -all: calculates both CU and RSCU for an input fasta or multifasta or a fetched ENST
+ -fetch: fetches a transcript (specifically the main CDS) from Ensembl, using the transcript identifier.
+ -fetchfile: same as fetch, but using a list of identifiers from a file, one per line.
+ -slice_domains: used to easily compare CU and RSCU between domains. It uses a file with "slice instructions" (find examples below).
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
To use the slice function (an example slice.csv is provided to show how the file needs to look like):

```bash
./coconut -cu -slice_domains  example.fasta slice.csv example.csv
```
To perform %MinMax calculations (no output name is required here):

```bash
./coconut -minmax Usage-num3.man multi_cds.fasta
```
And if you want to use a different window of codons:

```bash
./coconut -minmax  Usage-num3.man 10 multi_cds.fasta
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
| Fetch multiple versions of a transcript     | ⏰        |
| Slice domains from fetch                    | ⏰        |
| Slice multiple versions of a transcript & er| ⏰        |
| Fetch domains from InterPro or CATH         | ❌        |
| Incorrect length error                      | ❌        |
| Optional flags (no ATG, no STOP, etc)       | ❌        |
| MobiDB integration for disorder regions     | ❌        |

---



