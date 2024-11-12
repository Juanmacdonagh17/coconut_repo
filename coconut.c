

/////////////////////////
//        libs         //
/////////////////////////

// libraries that are needed
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h> 
#include <curl/curl.h> // guarda que acá en el gcc tengo que pasar manualmente donde esta (path de miniconda)

// path local para el gcc: -I/home/juan/miniconda3/include -L/home/juan/miniconda3/lib -lcurl


/////////////////////////
//  structures         //
/////////////////////////

// pointers
#define CODON_LENGTH 3
#define MAX_CODONS 64
#define MAX_LINE_LENGTH 1024 // maybe play around with this, so it can take longer sequences?
#define MAX_SEQUENCE_LENGTH 200000 
#define ID_LENGTH 512 // hopefully 512 is big enough, codonw tends to cut off the id
#define API_URL_FORMAT "https://rest.ensembl.org/sequence/id/%s?object_type=transcript;type=cds;content-type=text/x-fasta"  // dynamic ID URL
// here we could define a new API_URL_FORMAT where instead of one it gets all of the different versions of the transcirpt using multiple_sequences=1
// #define API_URL_FORMAT_MULT "https://rest.ensembl.org/sequence/id/%s?object_type=transcript;type=cds;content-type=text/x-fasta;multiple_sequences=1
// there needs to be a flag in main and then it should be treated as a multifasta within the fetch and multi fetch :)



// structures that we will use as output
typedef struct {
    char codon[CODON_LENGTH + 1];

    int count; // count
    int amino_acid_index; // for minmax

    float cu; // codon usage
    float rscu; // relative synonymous codon usage
    float usage; // for minmax
    float max; // for minmax
    float min; // for minmax
    float ave; // for minmax

} CodonCount;

// structure to handle HTTP response
struct MemoryStruct {
    char *memory;
    size_t size;
};

// structure for protein ID, codon count and nt seq
typedef struct {
    char id[ID_LENGTH]; 
    CodonCount counts[MAX_CODONS];
    char sequence[MAX_SEQUENCE_LENGTH]; // adjust size??
    } SequenceCodonCounts;

// structure for slicig and domain analysis
typedef struct {
    char transcript_id[ID_LENGTH];
    char domain_name[ID_LENGTH];
    int start;
    int end;
} SliceInstruction;

// minmax params
typedef struct {
    float usage;
    float max;
    float min;
    float ave;
    int amino_acid_index;
} CodonUsageStats;
CodonUsageStats codon_usage_stats[4][4][4]; 


// 64 codons
const char* codons[MAX_CODONS] = { 
    "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT",
    "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
    "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT",
    "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
    "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
    "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
    "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT",
    "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"
};

// Aa list, stops are coded as "*"
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

/////////////////////////
//     functions       //
/////////////////////////

////////////////////////
// request functions: //
////////////////////////

// using libcurl to write the fetched data into a memory buffer
static size_t WriteMemoryCallback(void *contents, size_t size, size_t nmemb, void *userp) {
    size_t realsize = size * nmemb;
    struct MemoryStruct *mem = (struct MemoryStruct *)userp;

    char *ptr = realloc(mem->memory, mem->size + realsize + 1);
    if(ptr == NULL) {
        // out of memory error
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
            fprintf(stderr, "curl_easy_perform() failed: %s\n", curl_easy_strerror(res)); // 17/10 every time i get this error?? 
                                                                                          // ENSEMBLE WTF I changed nothing and now worls? gotta keep looking
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

// retrives index of Aa of a codon
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

            // i think there is a bug here? specially with TTT. 

            // FIXED... there was an issue with how the table and the stops were ordered!! 
            
            return amino_acids[i] - 'A';  
            
        }
    }
    return -1; // error case if codon not found <- also this, I should try some wierd sequences.
               // also I think the handeling of the N's and the "weird" nt should be here ?  // 
               // at the countCodons functions the "N" are ignored
}

// converts codon to nt indices, for slicing :p 
void getCodonIndices(const char *codon, int *i, int *j, int *k) {
    char nucleotide_map[256] = {0};
    nucleotide_map['A'] = 0;
    nucleotide_map['C'] = 1;
    nucleotide_map['G'] = 2;
    nucleotide_map['T'] = 3;

    *i = nucleotide_map[(unsigned char)codon[0]];
    *j = nucleotide_map[(unsigned char)codon[1]];
    *k = nucleotide_map[(unsigned char)codon[2]];
}

////////////////////////
// numeric functions: //
////////////////////////

// starts counts of codons for each sequence
void initializeCodonCounts(SequenceCodonCounts *seq) {
    for (int i = 0; i < MAX_CODONS; i++) { // start all of the values at 0, maximum should be the length of the sequence of codons
        strcpy(seq->counts[i].codon, codons[i]);
        seq->counts[i].count = 0; // start count at 0
        seq->counts[i].cu = 0.0; // start CU at 0 :P 
        // maybe here too we should start RSCU at 0?
        seq->counts[i].rscu = 0.0; // yes? 
        
    }
    seq->sequence[0] = '\0'; // initialize sequence to empty string

}
// counts codons for each sequence
void countCodons(char *sequence, SequenceCodonCounts *seq, bool calculate_minmax) {
    int len = strlen(sequence);

    for (int i = 0; i <= len - CODON_LENGTH; i += CODON_LENGTH) {
        char codon[CODON_LENGTH + 1] = {0};
        strncpy(codon, &sequence[i], CODON_LENGTH);

        // check for 'N' in the codon
        if (strchr(codon, 'N') != NULL) {
            continue;
        }

        for (int j = 0; j < MAX_CODONS; j++) {
            if (strcmp(seq->counts[j].codon, codon) == 0) {
                seq->counts[j].count++;

                if (calculate_minmax) {
                    int idx_i, idx_j, idx_k;
                    getCodonIndices(codon, &idx_i, &idx_j, &idx_k);
                    if (idx_i >= 0 && idx_j >= 0 && idx_k >= 0) {
                        CodonUsageStats stats = codon_usage_stats[idx_i][idx_j][idx_k];
                        seq->counts[j].usage = stats.usage;
                        seq->counts[j].max = stats.max;
                        seq->counts[j].min = stats.min;
                        seq->counts[j].ave = stats.ave;
                        seq->counts[j].amino_acid_index = stats.amino_acid_index;
                    }
                }

                break;
            }
        }
    }
}
// calculates CU for each sequence
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
// calculates RSCU for each sequence
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
// MinMax functions:  //
////////////////////////

// reads codon usage data from a file
void readCodonUsageData(const char *filename) {
    FILE *infile = fopen(filename, "r");
    if (!infile) {
        fprintf(stderr, "Error opening codon usage file: %s\n", filename);
        exit(1);
    }

    int i, j, k, amino_index;
    float usage, max, min, ave;
    char line[MAX_LINE_LENGTH];

    while (fgets(line, sizeof(line), infile)) {
        if (sscanf(line, "%d %d %d %f %f %f %f %d", &i, &j, &k, &usage, &max, &min, &ave, &amino_index) == 8) {
            codon_usage_stats[i][j][k].usage = usage;
            codon_usage_stats[i][j][k].max = max;
            codon_usage_stats[i][j][k].min = min;
            codon_usage_stats[i][j][k].ave = ave;
            codon_usage_stats[i][j][k].amino_acid_index = amino_index;
        }
    }

    fclose(infile);
}

// calculates minmax
void calculateMinMax(SequenceCodonCounts *seq, int window_size, const char *output_filename) {
    FILE *minmax_output = fopen(output_filename, "a"); // open in append mode, so we can also proces mulit fastas now
    if (!minmax_output) {
        fprintf(stderr, "Error opening minmax output file: %s\n", output_filename);
        exit(1);
    }

    int len = strlen(seq->sequence);
    int num_codons = len / CODON_LENGTH;

    for (int j = 0; j <= num_codons - window_size; j++) {
        float sumusage = 0.0, summax = 0.0, summin = 0.0, sumave = 0.0;
        int valid_codons = 0;

        for (int k = 0; k < window_size; k++) {
            int idx = j + k;
            char codon[CODON_LENGTH + 1];
            strncpy(codon, &seq->sequence[idx * CODON_LENGTH], CODON_LENGTH);
            codon[CODON_LENGTH] = '\0';

            if (strchr(codon, 'N') != NULL) {
                continue;
            }

            int idx_i, idx_j, idx_k;
            getCodonIndices(codon, &idx_i, &idx_j, &idx_k);
            if (idx_i >= 0 && idx_j >= 0 && idx_k >= 0) {
                CodonUsageStats stats = codon_usage_stats[idx_i][idx_j][idx_k];
                if (stats.amino_acid_index < 1000) {
                    sumusage += stats.usage;
                    summax += stats.max;
                    summin += stats.min;
                    sumave += stats.ave;
                    valid_codons++;
                }
            }
        }

        if (valid_codons == 0) {
            continue;
        }

        sumusage /= valid_codons;
        summax /= valid_codons;
        summin /= valid_codons;
        sumave /= valid_codons;

        float minmax_percent = 0.0;
        if (sumusage > sumave && (summax - sumave) != 0) {
            minmax_percent = (sumusage - sumave) / (summax - sumave) * 100.0;
        } else if (sumusage < sumave && (sumave - summin) != 0) {
            minmax_percent = - (sumave - sumusage) / (sumave - summin) * 100.0;
        }

        // get codon and amino acid for current position
        char codon_out[CODON_LENGTH + 1];
        strncpy(codon_out, &seq->sequence[j * CODON_LENGTH], CODON_LENGTH);
        codon_out[CODON_LENGTH] = '\0';

        int amino_index = getAminoAcidIndex(codon_out); // INDEX FUNCTION IS NOT APPENDING THE CORRECT CODON!!! 
        char amino_name[4] = "Xaa";  // default unknown amino acid, maybe I can play a bit around here
        if (amino_index >= 0 && amino_index < MAX_CODONS) {
            if (amino_acids[amino_index] == '*') {
                strcpy(amino_name, "*aa");
            } else if (amino_index >= 0 && amino_index < 26) { // assuming 26 standard amino acids
                amino_name[0] = amino_acids[amino_index];
                // amino_name[1] = 'a';
                // amino_name[2] = 'a';
                amino_name[1] = '\0'; //3
            }
        }

        // include Sequence ID in the output
        fprintf(minmax_output, "%s,%d,%s,%s,%.2f,%.2f,%.2f,%.2f,%.2f\n",
                seq->id, j + 1, codon_out, amino_name,
                sumusage, summax, summin, sumave, minmax_percent);
    }

    fclose(minmax_output);
}


////////////////////////
// CSV functions:     //
////////////////////////


// write the header for the CSV
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

// write the rows for the CSV
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


////////////////////////
// slicer functions:  //
////////////////////////

// parse the file with the instructions for slicing the sequence
void parseSliceFile(const char *filename, SliceInstruction **instructions, int *count) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening slice file: %s\n", filename);
        exit(1);
    }

    char line[MAX_LINE_LENGTH];
    int instr_count = 0;
    int capacity = 100;
    SliceInstruction *temp_instructions = malloc(capacity * sizeof(SliceInstruction));
    if (!temp_instructions) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    while (fgets(line, sizeof(line), file)) {
        if (instr_count >= capacity) {
            capacity *= 2;
            temp_instructions = realloc(temp_instructions, capacity * sizeof(SliceInstruction));
            if (!temp_instructions) {
                fprintf(stderr, "Memory reallocation failed\n");
                exit(1);
            }
        }

        // remove newline characters
        line[strcspn(line, "\r\n")] = '\0';

        // split the line by commas (ID, domain_name, range)
        char *token = strtok(line, ",");
        if (token) strcpy(temp_instructions[instr_count].transcript_id, token);

        token = strtok(NULL, ",");
        if (token) strcpy(temp_instructions[instr_count].domain_name, token);

        token = strtok(NULL, ":");
        if (token) temp_instructions[instr_count].start = atoi(token) - 1;  // zero-index

        token = strtok(NULL, ":");
        if (token) temp_instructions[instr_count].end = atoi(token) - 1;    // zero-index

        // debugging: Print parsed instruction
        printf("Parsed Instruction %d: Transcript_ID=%s, Domain_Name=%s, Start=%d, End=%d\n",
               instr_count, temp_instructions[instr_count].transcript_id,
               temp_instructions[instr_count].domain_name,
               temp_instructions[instr_count].start, temp_instructions[instr_count].end);

        instr_count++;
    }

    fclose(file);
    *instructions = temp_instructions;
    *count = instr_count;
}



// be aware that you need to select the # of codon (that is also the corresponding position of the Aa in the sequence!)
// a seuqnece that has 30 nt, has 10 codons, and so 10 Aa. If you use "30" for the last Aa, it will raise an out of bounds error
// maybe add a flag for this?

// this part requeires some serious debuging 

void sliceAndWriteSequence(SequenceCodonCounts *seq, SliceInstruction *instructions, int instr_count, FILE *output, bool calculate_cu, bool calculate_rscu, bool silent, bool calculate_minmax, int window_size, const char *output_filename) {    for (int i = 0; i < instr_count; i++) {
        if (strcmp(seq->id, instructions[i].transcript_id) == 0) {

            printf("Matched Slice Instruction %d: Transcript_ID=%s\n", i, instructions[i].transcript_id);

            int start_pos = instructions[i].start * CODON_LENGTH;
            int end_pos = instructions[i].end * CODON_LENGTH;

            printf("Slicing: %s from %d to %d, positions %d to %d\n", instructions[i].transcript_id, instructions[i].start + 1, instructions[i].end + 1, start_pos, end_pos);

            // validate slice bounds
            if (start_pos < 0 || end_pos > (int)strlen(seq->sequence) || start_pos >= end_pos) {
                fprintf(stderr, "Invalid slicing bounds for %s: start=%d, end=%d\n", seq->id, instructions[i].start + 1, instructions[i].end + 1);
                continue;
            }

            int slice_length = end_pos - start_pos;
            if (slice_length <= 0) {
                fprintf(stderr, "Invalid slice length for %s\n", instructions[i].domain_name);
                continue;
            }

            // allocate memory for sliced_sequence dynamically
            char *sliced_sequence = malloc(slice_length + 1);
            if (!sliced_sequence) {
                fprintf(stderr, "Memory allocation failed for sliced_sequence\n");
                continue;
            }

            strncpy(sliced_sequence, seq->sequence + start_pos, slice_length);
            sliced_sequence[slice_length] = '\0';
            printf("Sliced sequence: %s\n", sliced_sequence);
            printf("\n");


            // process the sliced sequence
            SequenceCodonCounts slicedSeq;
            initializeCodonCounts(&slicedSeq);
            strncpy(slicedSeq.id, instructions[i].domain_name, sizeof(slicedSeq.id) - 1);
            slicedSeq.id[sizeof(slicedSeq.id) - 1] = '\0';
            strncpy(slicedSeq.sequence, sliced_sequence, sizeof(slicedSeq.sequence) - 1);
            slicedSeq.sequence[sizeof(slicedSeq.sequence) - 1] = '\0';

            countCodons(sliced_sequence, &slicedSeq, false); // last false is for the minmax 

            if (calculate_cu) {
                calculateCU(&slicedSeq);
            }
            if (calculate_rscu) {
                calculateRSCU(&slicedSeq);
            }

            writeCsvRow(output, &slicedSeq, calculate_cu, calculate_rscu);

            // free the allocated memory for sliced_sequence
            free(sliced_sequence);
            } else {
                printf("No match for Slice Instruction %d: Transcript_ID=%s with Sequence ID=%s\nRemember that ID's should match between files!!!!\n", i, instructions[i].transcript_id, seq->id);
        
        }
    }
}

////////////////////////
//   proccesing bulk  //
////////////////////////

void processSequence(const char *sequence_data, const char *sequence_id, FILE *output, bool calculate_cu, bool calculate_rscu, bool calculate_minmax, int window_size, const char *output_filename, bool silent, bool slice_sequence, const char *slice_filename) {
    // initialize sequence structure
    SequenceCodonCounts currentSequence;
    initializeCodonCounts(&currentSequence);
    strncpy(currentSequence.id, sequence_id, sizeof(currentSequence.id) - 1);
    currentSequence.id[sizeof(currentSequence.id) - 1] = '\0';

    // copy the sequence data
    strncpy(currentSequence.sequence, sequence_data, sizeof(currentSequence.sequence) - 1);
    currentSequence.sequence[sizeof(currentSequence.sequence) - 1] = '\0';

    // count codons and perform calculations
    countCodons(currentSequence.sequence, &currentSequence, calculate_minmax);

    if (calculate_cu) {
        calculateCU(&currentSequence);
    }
    if (calculate_rscu) {
        calculateRSCU(&currentSequence);
    }

    // write to CSV
    if (!silent) {
        if (calculate_cu || calculate_rscu) {
            writeCsvRow(stdout, &currentSequence, calculate_cu, calculate_rscu);
            printf("\n");
        }   
        
    }
    
    if (calculate_cu || calculate_rscu) {
       // writeCsvHeader(output, calculate_cu, calculate_rscu);
        writeCsvRow(output, &currentSequence, calculate_cu, calculate_rscu); // the thing with this is that you also don't get the codon count... ehhh do we want it??
    }   

    // perform min-max calculations
    if (calculate_minmax) {
       // char minmax_output_filename[1024];
        //printf("Minmax output filename: %s\n", minmax_output_filename);

        // HERE there are some issues with how the output file is named. 
        // the output gets spaces and a binch of thins, it should be better handeled. 
        // 
       // snprintf(minmax_output_filename, sizeof(minmax_output_filename), "%s_minmax.out", currentSequence.id);
       // calculateMinMax(&currentSequence, window_size, minmax_output_filename);
        calculateMinMax(&currentSequence, window_size, output_filename);

    }

    // slicing functionality
    if (slice_sequence) {
        // parse the slice file if not already parsed
        static SliceInstruction *instructions = NULL;
        static int slice_count = 0;
        static bool slices_parsed = false;
        if (!slices_parsed) {
            parseSliceFile(slice_filename, &instructions, &slice_count);
            slices_parsed = true;
        }

        // process each sliced sequence
        sliceAndWriteSequence(&currentSequence, instructions, slice_count, output, calculate_cu, calculate_rscu, silent, calculate_minmax, window_size, output_filename);
    }
}


// bool isEnsemblServerOnline() {
//     CURL *curl;
//     CURLcode res;
//     long response_code = 0;
//     bool server_online = false;

//     curl = curl_easy_init();
//     if(curl) {
//         // Set the URL to the Ensembl REST API base URL
//         curl_easy_setopt(curl, CURLOPT_URL, "https://rest.ensembl.org/");
//         // Perform a HEAD request
//         curl_easy_setopt(curl, CURLOPT_NOBODY, 1L);
//         // Set a timeout for the request
//         curl_easy_setopt(curl, CURLOPT_TIMEOUT, 10L);

//         res = curl_easy_perform(curl);

//         if(res == CURLE_OK) {
//             // Check the HTTP response code
//             curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &response_code);
//             if(response_code == 200) {
//                 server_online = true;
//             }
//         } else {
//             fprintf(stderr, "Failed to connect to Ensembl server: %s\n", curl_easy_strerror(res));
//         }

//         curl_easy_cleanup(curl);
//     }
//     return server_online;
// }


////////////////////////
//       main         //
////////////////////////

int main(int argc, char *argv[]) {
    bool calculate_cu = false;
    bool calculate_rscu = false;
    bool silent = false;  // silent flag so you don't get all the printf or fprintf:
    bool fetch_sequence = false;
    bool fetch_from_file = false; // fetch from a .txt
    bool slice_sequence = false;
    bool calculate_minmax = false;

    int window_size = 18;
    int slice_count = 0;

    char *codon_usage_filename = NULL;
    char *slice_filename = NULL;  // path to slice file
    char *input_filename = NULL;
    char *output_filename = NULL;
    char *protein_id = NULL;
    char *fetchfile = NULL; // fetch from a .txt

    // variable declarations for slicing
    SliceInstruction *instructions = NULL;


    for (int i = 1; i < argc; i++) {

        if (strcmp(argv[i], "-cu") == 0) {
            calculate_cu = true;

        } else if (strcmp(argv[i], "-rscu") == 0) {
            calculate_rscu = true;

        } else if (strcmp(argv[i], "-all") == 0) {
            calculate_cu = true;
            calculate_rscu = true;
        } else if (strcmp(argv[i], "-silent") == 0) {
            silent = true;  // silent mode

        } else if (strcmp(argv[i], "-fetch") == 0) {
            if (i + 2 < argc) {
                fetch_sequence = true;
                protein_id = argv[i + 1];
                output_filename = argv[i + 2];
                i += 2;  } // adjust 'i' after consuming arguments 

        } else if (strcmp(argv[i], "-fetchfile") == 0) {
            if (i + 2 < argc) {
                fetch_from_file = true;
                fetchfile = argv[i + 1];
                output_filename = argv[i + 2];
                i += 2; }

        } else if (strcmp(argv[i], "-slice_domains") == 0) {
            if (i + 1 < argc) {
                slice_sequence = true;
                slice_filename = argv[i + 1];
                i += 1; } // get the path to the slice file

        } else if (strcmp(argv[i], "-minmax") == 0) {
            if (i + 2 < argc) {
                calculate_minmax = true;
                codon_usage_filename = argv[++i];   // table file
                input_filename = argv[++i];         // fasta file

                // check if the next argument is a number (window size)
                if (i + 1 < argc && isdigit(argv[i + 1][0])) {
                    window_size = atoi(argv[++i]);
                } else {

                 fprintf(stderr, "Using default codon window size: 18\n");
                }

            } else {
                fprintf(stderr, "Error: Codon usage file and input fasta file must be specified with -minmax flag.\n");
                return 1;
            }

        
        } else if (strcmp(argv[i], "-help") == 0) {

            printf("Usage: %s [options] <input.fasta> [output.csv]\n", argv[0]);
            printf("Options:\n");
            printf("  -cu\tCalculate Codon Usage (CU)\n");
            printf("  -rscu\tCalculate Relative Synonymous Codon Usage (RSCU)\n");
            printf("  -all\tCalculate CU and RSCU\n");
            printf("  -silent\tHide outputs from the console\n");
            printf("  -fetch\tFetch sequence from Ensembl using transcript ID (ENST) (requires protein_id and output_filename)\n");
            printf("  -fetchfile\tFetch multiple sequences from a text file (requires fetchfile and output_filename)\n");
            printf("  -slice_domains\tSlice fasta into domains using a CSV-style file (requires slice_filename)\n");
            printf("  -minmax <codon_usage_file> [window_size]\tCalculate min-max percentage over specified window size (default 18).\n");
            printf("\n");
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

    ///////////////////////////////////////////
    // errors for arguments or files missing //
    ///////////////////////////////////////////
    if (!input_filename) {
        fprintf(stderr, "Input filename is required.\n");
        return 1;
    }
    // if ((!input_filename || !output_filename) && !fetch_sequence && !fetch_from_file) {
    //     fprintf(stderr, "Unkwon command:\n Usage: %s [-cu] [-rscu] [-all] [-silent] [-fetch] [-fetchfile] [-slice_domains] [-minmax <codon_usage_file> [window_size]] [-help] <slice.csv> <input.fasta> <output.csv>\n", argv[0]);
    //     return 1;
    // }
    // check for output file in both fetch and normal input cases
    if ((calculate_cu || calculate_rscu || slice_sequence || fetch_sequence || fetch_from_file) && !output_filename) {
        fprintf(stderr, "Output filename is required for the selected options.\n");
        return 1;
    }
    // if (!output_filename) {
    //     fprintf(stderr, "Output filename is required.\n");
    //     return 1;
    // }
 
    FILE *output = NULL;
    if (calculate_cu || calculate_rscu || slice_sequence || fetch_sequence || fetch_from_file) {
        output = fopen(output_filename, "w");
        if (!output) {
            if (!silent) {
                fprintf(stderr, "Error opening output file: %s\n", output_filename);
            }
            return 1;
        }

        // write CSV header for CU and RSCU
        if (calculate_cu || calculate_rscu) {
            writeCsvHeader(output, calculate_cu, calculate_rscu);
        }
    }
    if (calculate_minmax) {
        if (!codon_usage_filename) {
            fprintf(stderr, "Error: Codon usage file must be specified with -minmax flag.\n");
            fclose(output);
            return 1;
        }
        // Read codon usage data
        readCodonUsageData(codon_usage_filename);
        printf("Codon usage data loaded successfully.\n");

        // Generate output filename based on input_filename
        char *dot = strrchr(input_filename, '.');
        char *base_filename = NULL;
        if (dot) {
            base_filename = strndup(input_filename, dot - input_filename);
        } else {
            base_filename = strdup(input_filename);
        }
        output_filename = malloc(strlen(base_filename) + 13); // '_minmax.out' + null terminator
        sprintf(output_filename, "%s_minmax.out", base_filename);
        free(base_filename);

        // Open the output file to write the header
        FILE *minmax_output = fopen(output_filename, "w");
        if (!minmax_output) {
            fprintf(stderr, "Error opening minmax output file: %s\n", output_filename);
            fclose(output);
            return 1;
        }
        // Write header
        fprintf(minmax_output, "ID,position,codon,aminoacid,usage,max,min,ave,minmax\n");
        fclose(minmax_output);
    }   
 

    char *sequence_data = NULL;




    // args //



    if (fetch_sequence) { // as of 17/10 i can't use the fetch option wtf did I broke NOW PLEASE 

        // if (!isEnsemblServerOnline()) {
        //     fprintf(stderr, "Ensembl server is not available. Please try again later.\n");
        //     fclose(output);
        //     return 1;
        // }
        // fetch a single sequence
        if (!silent) {
            printf("Fetching sequence for protein ID: %s\n", protein_id);
        }

        char *sequence_data = fetchProteinSequence(protein_id);
        if (!sequence_data) {
            fprintf(stderr, "Failed to fetch sequence for protein ID: %s\n", protein_id);
            fclose(output);
            return 1;
        }

        if (!silent) {
            printf("Fetched Sequence:\n%s\n", sequence_data);
        }

        // process the fetched sequence
        processSequence(sequence_data, protein_id, output, calculate_cu, calculate_rscu, calculate_minmax, window_size, output_filename, silent, slice_sequence, slice_filename);

        free(sequence_data);

    } else if (fetch_from_file) {
        // fetch multiple sequences from a file
        FILE *file = fopen(fetchfile, "r");
        if (!file) {
            fprintf(stderr, "Error opening file: %s\n", fetchfile);
            fclose(output);
            return 1;
        }

        char enst_id[MAX_LINE_LENGTH];
        while (fgets(enst_id, sizeof(enst_id), file)) {
            enst_id[strcspn(enst_id, "\r\n")] = '\0';  // remove newline characters

            if (!silent) {
                printf("Fetching sequence for protein ID: %s\n", enst_id);
            }

            char *sequence_data = fetchProteinSequence(enst_id);
            if (!sequence_data) {
                if (!silent) {
                    fprintf(stderr, "Failed to fetch sequence for ENST ID: %s. Skipping...\n", enst_id);
                }
                continue;
            }

            if (!silent) {
                printf("Fetched Sequence:\n%s\n", sequence_data);
            }

            // process the fetched sequence
            processSequence(sequence_data, enst_id, output, calculate_cu, calculate_rscu, calculate_minmax, window_size, output_filename, silent, slice_sequence, slice_filename);

            free(sequence_data);
        }
        fclose(file);

    } else if (input_filename) {
        // read sequences from a FASTA file
        FILE *input = fopen(input_filename, "r");
        if (!input) {
            if (!silent) {
                fprintf(stderr, "Failed to open input file: %s\n", input_filename);
            }
            fclose(output);
            return 1;
        }

        SequenceCodonCounts currentSequence;
        bool active = false;
        char line[MAX_LINE_LENGTH];

        while (fgets(line, sizeof(line), input)) {
            line[strcspn(line, "\r\n")] = '\0';  // remove newline

            if (line[0] == '>') {  // header line
                if (active) {
                    // process the previous sequence
                    processSequence(currentSequence.sequence, currentSequence.id, output, calculate_cu, calculate_rscu, calculate_minmax, window_size, output_filename, silent, slice_sequence, slice_filename);
                }
                // start a new sequence
                strncpy(currentSequence.id, line + 1, sizeof(currentSequence.id) - 1);
                currentSequence.id[sizeof(currentSequence.id) - 1] = '\0';
                initializeCodonCounts(&currentSequence);
                active = true;
            } else if (active && strlen(line) > 0) {
                // accumulate sequence data
                for (int i = 0; line[i]; i++) {
                    line[i] = toupper(line[i]);  // convert to uppercase
                }
                strcat(currentSequence.sequence, line);
            }
        }

        // process the last sequence
        if (active) {
            processSequence(currentSequence.sequence, currentSequence.id, output, calculate_cu, calculate_rscu, calculate_minmax, window_size, output_filename, silent, slice_sequence, slice_filename);
        }

        fclose(input);
    }

    // free and close stuff 

    if (output) {
        fclose(output);
    }

    if (output_filename) free(output_filename);
    return 0;

}







