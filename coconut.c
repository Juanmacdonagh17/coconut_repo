#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h> 
#include <curl/curl.h> // guarda que acá en el gcc tengo que pasar manualmente donde esta (path de miniconda)


// this are the pointers, i have no idea about pointers
#define CODON_LENGTH 3
#define MAX_CODONS 64
#define MAX_LINE_LENGTH 1024 // maybe play around with this, so it can take longer sequences?
#define ID_LENGTH 512 // hopefully 512 is big enough, codonw tends to cut off the id
#define API_URL_FORMAT "https://rest.ensembl.org/sequence/id/%s?object_type=transcript;type=cds;content-type=text/x-fasta"  // dynamic ID URL

// numbers that we will use as output
typedef struct {
    char codon[CODON_LENGTH + 1];
    int count; // count
    float cu; // codon usage
    float rscu; // relative synonymous codon usage
} CodonCount;

// structure to handle HTTP response
struct MemoryStruct {
    char *memory;
    size_t size;
};

typedef struct {
    char id[ID_LENGTH]; 
    CodonCount counts[MAX_CODONS];
} SequenceCodonCounts;

const char* codons[MAX_CODONS] = { // 64 codons, nothing wierd here 
    "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT",
    "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
    "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT",
    "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
    "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
    "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
    "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT",
    "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"
};

// correct order of Aa, stops are coded as "*":
const char amino_acids[MAX_CODONS] = {
    'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T',
    'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S',
    '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'
};

////////////////////////
// request functions: //
////////////////////////

static size_t WriteMemoryCallback(void *contents, size_t size, size_t nmemb, void *userp) {
    size_t realsize = size * nmemb;
    struct MemoryStruct *mem = (struct MemoryStruct *)userp;

    char *ptr = realloc(mem->memory, mem->size + realsize + 1);
    if(ptr == NULL) {
        // Out of memory error
        return 0;
    }

    mem->memory = ptr;
    memcpy(&(mem->memory[mem->size]), contents, realsize);
    mem->size += realsize;
    mem->memory[mem->size] = 0;

    return realsize;
}

// fetch protein sequence from Ensembl
char* fetchProteinSequence(const char* protein_id) {
    CURL *curl_handle;
    CURLcode res;

    struct MemoryStruct chunk;
    chunk.memory = malloc(1);  // initialize memory
    chunk.size = 0;

    // dynamic URL using the protein ENST id
    char url[256];
    snprintf(url, sizeof(url), API_URL_FORMAT, protein_id);

    curl_global_init(CURL_GLOBAL_ALL);  

    // libcurl -> be aware that it has to be installed prior (it can be also present in conda or miniconda)
    curl_handle = curl_easy_init();

    if(curl_handle) {
        curl_easy_setopt(curl_handle, CURLOPT_URL, url);
        curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, WriteMemoryCallback);
        curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, (void *)&chunk);
        curl_easy_setopt(curl_handle, CURLOPT_USERAGENT, "libcurl-agent/1.0");
       
        res = curl_easy_perform(curl_handle);
       
        if(res != CURLE_OK) {
            fprintf(stderr, "curl_easy_perform() failed: %s\n", curl_easy_strerror(res));
            return NULL;
        }

        curl_easy_cleanup(curl_handle);
    }

    curl_global_cleanup();

    return chunk.memory;  // return fetched sequence
}

////////////////////////
// index functions:   //
////////////////////////

int getAminoAcidIndex(char* codon) {
    for (int i = 0; i < MAX_CODONS; i++) { 

        // check if the current codon in the array matches the input codon using string comparison

        if (strcmp(codons[i], codon) == 0) {

            if (amino_acids[i] == '*') {
                return 25; // use 25 as a special index for stop codons
            }
            
            // if a match is found, compute the index:
            // subtract 'A' from the corresponding amino acid character to get an index.
            // this index is calculated based on the position of the amino acid character in the alphabet,
            // where 'A' is 0, 'B' is 1, ..., 'Z' is 25.

            // I think there is a bug here? specially with TTT. 

            // FIXED... there was an issue with how the table and the stops were ordered!! 
            
            return amino_acids[i] - 'A';  
            
        }
    }
    return -1; // error case if codon not found <- also this, I should try some wierd sequences.
               // also I think the handeling of the N's and the "weird" nt should be here ?  // 
               // at the countCodons functions the "N" are ignored
}

////////////////////////
// numeric functions: //
////////////////////////


void initializeCodonCounts(SequenceCodonCounts *seq) {
    for (int i = 0; i < MAX_CODONS; i++) { // start all of the values at 0, maximum should be the length of the sequence of codons
        strcpy(seq->counts[i].codon, codons[i]);
        seq->counts[i].count = 0; // start count at 0
        seq->counts[i].cu = 0.0; // start CU at 0 :P 
        // maybe here too we should start RSCU at 0?
    }
}

void countCodons(char *sequence, SequenceCodonCounts *seq) {
   int len = strlen(sequence);
   // int completeCodonLength = len - (len % CODON_LENGTH); // Ensure only complete codons are processed

    for (int i = 0; i <= len - CODON_LENGTH; i += CODON_LENGTH) {
        char codon[CODON_LENGTH + 1] = {0};
        strncpy(codon, &sequence[i], CODON_LENGTH);

        // check if the codon contains 'N'
        if (strchr(codon, 'N') != NULL) {
            continue; // skip this codon if it contains 'N' <- maybe it will be usefull to keep track of this?
        }

        for (int j = 0; j < MAX_CODONS; j++) {
            if (strcmp(seq->counts[j].codon, codon) == 0) {
                seq->counts[j].count++;
                break;
            }
        }
    }
}

void calculateCU(SequenceCodonCounts *seq) {
    int amino_totals[26] = {0}; // array to store total counts for each amino acid

    for (int i = 0; i < MAX_CODONS; i++) {
        int amino_index = getAminoAcidIndex(seq->counts[i].codon);
        if (amino_index >= 0) { // ensure the index is valid
            amino_totals[amino_index] += seq->counts[i].count;
        }
    }

    for (int i = 0; i < MAX_CODONS; i++) {
        int amino_index = getAminoAcidIndex(seq->counts[i].codon);
        if (amino_index >= 0 && amino_totals[amino_index] > 0) { // avoid div by 0
            seq->counts[i].cu = (float)seq->counts[i].count / amino_totals[amino_index];
        }
    }
}

void calculateRSCU(SequenceCodonCounts *seq) {
    int amino_totals[26] = {0}; // 26 letters in the alphabet, assuming one amino acid per letter (with stops)
    int synonymous_counts[MAX_CODONS] = {0}; // array to store counts of synonymous codons

    for (int i = 0; i < MAX_CODONS; i++) {
        int amino_index = getAminoAcidIndex(seq->counts[i].codon);
        if (amino_index >= 0) { // ensure the index is valid
            amino_totals[amino_index] += seq->counts[i].count;
            synonymous_counts[amino_index]++;
        }
    }

    for (int i = 0; i < MAX_CODONS; i++) {
        int amino_index = getAminoAcidIndex(seq->counts[i].codon);
        if (amino_index >= 0 && amino_totals[amino_index] > 0) { // avoid division by zero
            seq->counts[i].rscu = (float)seq->counts[i].count / (amino_totals[amino_index] / (float)synonymous_counts[amino_index]);
        } else {
            seq->counts[i].rscu = 0.0; // ensure RSCU is zero if totals are zero
        }
    }
}


////////////////////////
// CSV functions:     //
////////////////////////

void writeCsvHeader(FILE *file, bool includeCU, bool includeRSCU) {
    fprintf(file, "ID");
    for (int i = 0; i < MAX_CODONS; i++) {
        fprintf(file, ",%s", codons[i]);
        if (includeCU) {
            fprintf(file, ",%s_CU", codons[i]);
        }
        if (includeRSCU) {
            fprintf(file, ",%s_RSCU", codons[i]);
        }        
    }
    fprintf(file, "\n");
}

void writeCsvRow(FILE *file, SequenceCodonCounts *seq, bool includeCU, bool includeRSCU) {
    fprintf(file, "%s", seq->id);
    for (int i = 0; i < MAX_CODONS; i++) {
        fprintf(file, ",%d", seq->counts[i].count);
        if (includeCU) {
            fprintf(file, ",%.3f", seq->counts[i].cu);
        }
        if (includeRSCU) {
            fprintf(file, ",%.3f", seq->counts[i].rscu);
        }
    }
    fprintf(file, "\n");
}


// flag handeling, file name, etc 

int main(int argc, char *argv[]) {
    bool calculate_cu = false;
    bool calculate_rscu = false;
    bool silent = false;  // silent flag so you don't get all the printf or fprintf:
    bool fetch_sequence = false;
    bool fetch_from_file = false; // fetch from a .txt

    char *input_filename = NULL;
    char *output_filename = NULL;
    char *protein_id = NULL;
    char *fetchfile = NULL; // fetch from a .txt


    for (int i = 1; i < argc; i++) {

        if (strcmp(argv[i], "-cu") == 0) {
            calculate_cu = true;

        } else if (strcmp(argv[i], "-rscu") == 0) {
            calculate_rscu = true;

        } else if (strcmp(argv[i], "-silent") == 0) {
            silent = true;  // silent mode

        } else if (strcmp(argv[i], "-fetch") == 0 && i + 1 < argc) {
            fetch_sequence = true;
            protein_id = argv[++i];
            output_filename = argv[++i];  // make sure to capture the output filename from the fetch

        } else if (strcmp(argv[i], "-fetchfile") == 0 && i + 1 < argc) {
            fetch_from_file = true;
            fetchfile = argv[++i];
            output_filename = argv[++i];

        } else if (strcmp(argv[i], "-all") == 0) {
            calculate_cu = true;
            calculate_rscu = true;

        } else if (strcmp(argv[i], "-help") == 0) {

            printf("Usage: %s [options] <input.fasta> <output.csv>\n", argv[0]);
            printf("Options:\n");
            printf("  -cu\tCalculate Codon Usage (CU)\n");
            printf("  -rscu\tCalculate Relative Synonymous Codon Usage (RSCU)\n");
            printf("  -all\tCalculate CU and RSCU\n");
            printf("  -silent\tHide outputs from the console\n");
            printf("  -fetch\tFetch sequence from Ensembl using transcript ID (ENST) (no <input.fasta>)\n");
            printf("  -fetchfile\tFetch multiple sequences from a text file (ENSTs)\n");

            printf("\n");
            printf("\n");
            printf("Developed by Juan Mac Donagh at JGU and UNQ\n");
            printf("Contact: macjuan17@gmail.com\n");

            return 0;

        } else if (!input_filename) {
            input_filename = argv[i];
        } else if (!output_filename) {
            output_filename = argv[i];
        }
    }

    if ((!input_filename || !output_filename) && !fetch_sequence && !fetch_from_file) {
        fprintf(stderr, "Usage: %s [-cu] [-rscu] [-all] [-silent] [-fetch] [-fetchfile] [-help] <input.fasta> <output.csv>\n", argv[0]);
        return 1;
    }

    // check for output file in both fetch and normal input cases
    if (!output_filename) {
        fprintf(stderr, "Output filename is required.\n");
        return 1;
    }

    FILE *output = fopen(output_filename, "w");
    if (!output) {
        if (!silent) {
            fprintf(stderr, "Error opening input file: %s\n", input_filename);
        }
        return 1;
    }

    // write CSV header -> there is an issue here with the version of the ENST when doing a fetch i believe
    writeCsvHeader(output, calculate_cu, calculate_rscu);

    char *sequence_data = NULL;

    if (fetch_sequence) {

        if (!silent) {
            printf("Fetching sequence for protein ID: %s\n", protein_id);
        }

        sequence_data = fetchProteinSequence(protein_id);
        if (!sequence_data) {
            fprintf(stderr, "Failed to fetch sequence for protein ID: %s\n", protein_id);
            fclose(output);  // close the file before returning
            return 1;
        }

        if (!silent) {
            printf("Fetched Sequence:\n%s\n", sequence_data);
        }

        // process the fetched sequence
        SequenceCodonCounts currentSequence;
        initializeCodonCounts(&currentSequence);
        strncpy(currentSequence.id, protein_id, sizeof(currentSequence.id) - 1);
        currentSequence.id[sizeof(currentSequence.id) - 1] = '\0';

        countCodons(sequence_data, &currentSequence);

        if (calculate_cu) {
            calculateCU(&currentSequence);
        }
        if (calculate_rscu) {
            calculateRSCU(&currentSequence);
        }

        // write fetched sequence to CSV file
        if (!silent) {
            writeCsvRow(stdout, &currentSequence, calculate_cu, calculate_rscu);  // Print to console
        }
        writeCsvRow(output, &currentSequence, calculate_cu, calculate_rscu);
        
        free(sequence_data);  // free the fetched sequence memory
    


    } else if (fetch_from_file) {  // NEW: fetch multiple sequences from file
    FILE *file = fopen(fetchfile, "r");
    if (!file) {
        fprintf(stderr, "Error opening file: %s\n", fetchfile);
        fclose(output);
        return 1;
    }

    char enst_id[MAX_LINE_LENGTH];
    while (fgets(enst_id, sizeof(enst_id), file)) {
        enst_id[strcspn(enst_id, "\r\n")] = '\0';  // Remove newline characters

        char *sequence_data = fetchProteinSequence(enst_id);

        printf("Fetching sequence for protein ID: %s\n", enst_id);
        if (!silent) {
            printf("Fetched Sequence:\n%s\n", sequence_data);
        }
        if (!sequence_data) {
            if (!silent) fprintf(stderr, "Failed to fetch sequence for ENST ID: %s\n", enst_id);
            continue;  // Skip this ENST and proceed to the next
        }

        
        SequenceCodonCounts currentSequence;
        initializeCodonCounts(&currentSequence);
        strncpy(currentSequence.id, enst_id, sizeof(currentSequence.id) - 1);
        countCodons(sequence_data, &currentSequence);

        if (calculate_cu) {
            calculateCU(&currentSequence);
        }
        if (calculate_rscu) {
            calculateRSCU(&currentSequence);
        }
        if (!silent) {
        writeCsvRow(stdout, &currentSequence, calculate_cu, calculate_rscu);
        }
        writeCsvRow(output, &currentSequence, calculate_cu, calculate_rscu);

        free(sequence_data);
    }
    fclose(file);


    
    } else if (input_filename) {
        FILE *input = fopen(input_filename, "r");
        if (!input) {
            if (!silent) {
                fprintf(stderr, "Failed to open sequence for protein ID: %s\n", protein_id); // esto esta mal? ya estamos en el input
            }
            fclose(output);  // close the file before returning
            return 1;
        }

        // process the FASTA file (input, not fetched)
        SequenceCodonCounts currentSequence;
        bool active = false;
        char line[MAX_LINE_LENGTH];
        
        while (fgets(line, sizeof(line), input)) {
            line[strcspn(line, "\r\n")] = '\0'; // remove newline

            if (line[0] == '>') { // header line
                if (active) {
                    if (calculate_cu) {
                        calculateCU(&currentSequence);
                    }
                    if (calculate_rscu) {
                        calculateRSCU(&currentSequence);
                    }
                    // write processed sequence to CSV file
                    if (!silent) {
                        fprintf(stderr, "Procesing protein ID: %s:\n", currentSequence.id);
                        writeCsvRow(stdout, &currentSequence, calculate_cu, calculate_rscu);  // Print to console
                        fprintf(stderr, "\n");
                    }
                    writeCsvRow(output, &currentSequence, calculate_cu, calculate_rscu);  // write to CSV file
                }
                strncpy(currentSequence.id, line + 1, sizeof(currentSequence.id) - 1); // copy ID, trimming '>'
                currentSequence.id[strcspn(currentSequence.id, "\n")] = '\0'; // remove newline
                initializeCodonCounts(&currentSequence);
                active = true;
            } else if (active && strlen(line) > 0) {
                for (int i = 0; line[i]; i++) {
                    line[i] = toupper(line[i]); // convert to uppercase
                }
                countCodons(line, &currentSequence);
            }
        }

        if (active) { // final sequence processing
            if (calculate_cu) {
                calculateCU(&currentSequence);
            }
            if (calculate_rscu) {
                calculateRSCU(&currentSequence);
            }
        if (!silent) { // this only makes sure that it prints the ID and the process for the last fasta in a multi fasta. it seats here, out of the while loop
            fprintf(stderr, "Procesing protein ID: %s:\n", currentSequence.id);
            writeCsvRow(stdout, &currentSequence, calculate_cu, calculate_rscu);  // print to console
            fprintf(stderr, "\n");
        }
        writeCsvRow(output, &currentSequence, calculate_cu, calculate_rscu);  // write to CSV file
        }

        fclose(input);
    }

    fclose(output);  // close the output file after writing

    return 0;
}







