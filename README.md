***Coconut is currently being developed.
Be aware of bugs, especially tyranids.***

---

![alt text](https://github.com/Juanmacdonagh17/coconut_repo/blob/main/logo/Natural-fresh-coconut-illustration-Premium-Vector-PNG.png)

# Coconut: 
## A suite for transcripts and codon calculations
---

### Instalation 

To install, you need to have the curl.h library.
The simplest way is installing Miniconda3 (https://docs.anaconda.com/miniconda/).

After that, just run the following command:

```bash
gcc -o coconut coconut.c -I/home/path/to/miniconda3/include -L/home/path/to/miniconda3/lib -lcurl
```
Replace /path/to with the folder where Conda was installed.

To run coconut.c, just execute ./coconut in the command line.
Using the -help flag will provide a list of all available commands.

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
| Fetch multiple versions of a transcript     | ❌        |
| Slice domains from fetch                    | ❌        |
| Slice multiple versions of a transcript & er| ❌        |
| Calculate MinMax                            | ❌        |
| Fetch domains from InterPro                 | ❌        |
| Incorrect length error                      | ❌        |
| Optional flags (no ATG, no STOP, etc        | ❌        |
| Docker image                                | ❌        |




