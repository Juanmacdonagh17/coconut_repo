#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


#define NUM_CODONS 64

typedef struct {
    char codon[4];      
    char aa_name[4];     
    int  aa_number;      
    int  base_code[3];   
    long count;          
    double frequency;    
} CodonInfo;

/* 
 *   A=0, C=1, G=2, T=3
 */
int base_to_code(char b) {
    switch (toupper((unsigned char)b)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        //default:  return -1; // for ambiguous or invalid
    }
}


CodonInfo codon_table[NUM_CODONS] = {
    // Codon,  AA,  AA#, [base_code0, base_code1, base_code2]
    {"GCT", "Ala", 1, {2,1,3}, 0, 0.0},
    {"GCC", "Ala", 1, {2,1,1}, 0, 0.0},
    {"GCA", "Ala", 1, {2,1,0}, 0, 0.0},
    {"GCG", "Ala", 1, {2,1,2}, 0, 0.0},

    {"CGT", "Arg", 2, {1,2,3}, 0, 0.0},
    {"CGC", "Arg", 2, {1,2,1}, 0, 0.0},
    {"CGA", "Arg", 2, {1,2,0}, 0, 0.0},
    {"CGG", "Arg", 2, {1,2,2}, 0, 0.0},
    {"AGA", "Arg", 2, {0,2,0}, 0, 0.0},
    {"AGG", "Arg", 2, {0,2,2}, 0, 0.0},

    {"AAT", "Asn", 3, {0,0,3}, 0, 0.0},
    {"AAC", "Asn", 3, {0,0,1}, 0, 0.0},

    {"GAT", "Asp", 4, {2,0,3}, 0, 0.0},
    {"GAC", "Asp", 4, {2,0,1}, 0, 0.0},

    {"TGT", "Cys", 5, {3,2,3}, 0, 0.0},
    {"TGC", "Cys", 5, {3,2,1}, 0, 0.0},

    {"CAA", "Gln", 6, {1,0,0}, 0, 0.0},
    {"CAG", "Gln", 6, {1,0,2}, 0, 0.0},

    {"GAA", "Glu", 7, {2,0,0}, 0, 0.0},
    {"GAG", "Glu", 7, {2,0,2}, 0, 0.0}, 

    {"GGA", "Gly", 8, {2,2,0}, 0, 0.0},
    {"GGC", "Gly", 8, {2,2,1}, 0, 0.0},
    {"GGG", "Gly", 8, {2,2,2}, 0, 0.0},
    {"GGT", "Gly", 8, {2,2,3}, 0, 0.0},

    {"CAT", "His", 9, {1,0,3}, 0, 0.0},
    {"CAC", "His", 9, {1,0,1}, 0, 0.0},

    {"ATT", "Ile", 10, {0,3,3}, 0, 0.0},
    {"ATC", "Ile", 10, {0,3,1}, 0, 0.0},
    {"ATA", "Ile", 10, {0,3,0}, 0, 0.0},

    {"TTA", "Leu", 11, {3,3,0}, 0, 0.0},
    {"TTG", "Leu", 11, {3,3,2}, 0, 0.0},
    {"CTT", "Leu", 11, {1,3,3}, 0, 0.0},
    {"CTC", "Leu", 11, {1,3,1}, 0, 0.0},
    {"CTA", "Leu", 11, {1,3,0}, 0, 0.0},
    {"CTG", "Leu", 11, {1,3,2}, 0, 0.0},

    {"AAA", "Lys", 12, {0,0,0}, 0, 0.0},
    {"AAG", "Lys", 12, {0,0,2}, 0, 0.0},

    {"ATG", "Met", 13, {0,3,2}, 0, 0.0},

    {"TTT", "Phe", 14, {3,3,3}, 0, 0.0},
    {"TTC", "Phe", 14, {3,3,1}, 0, 0.0},

    {"CCT", "Pro", 15, {1,1,3}, 0, 0.0},
    {"CCC", "Pro", 15, {1,1,1}, 0, 0.0},
    {"CCA", "Pro", 15, {1,1,0}, 0, 0.0},
    {"CCG", "Pro", 15, {1,1,2}, 0, 0.0},

    {"TCT", "Ser", 16, {3,1,3}, 0, 0.0},
    {"TCC", "Ser", 16, {3,1,1}, 0, 0.0},
    {"TCA", "Ser", 16, {3,1,0}, 0, 0.0},
    {"TCG", "Ser", 16, {3,1,2}, 0, 0.0},
    {"AGT", "Ser", 16, {0,2,3}, 0, 0.0},
    {"AGC", "Ser", 16, {0,2,1}, 0, 0.0},

    {"ACT", "Thr", 17, {0,1,3}, 0, 0.0},
    {"ACC", "Thr", 17, {0,1,1}, 0, 0.0},
    {"ACA", "Thr", 17, {0,1,0}, 0, 0.0},
    {"ACG", "Thr", 17, {0,1,2}, 0, 0.0},

    {"TGG", "Trp", 18, {3,2,2}, 0, 0.0},

    {"TAT", "Tyr", 19, {3,0,3}, 0, 0.0},
    {"TAC", "Tyr", 19, {3,0,1}, 0, 0.0},

    {"GTT", "Val", 20, {2,3,3}, 0, 0.0},
    {"GTC", "Val", 20, {2,3,1}, 0, 0.0},
    {"GTA", "Val", 20, {2,3,0}, 0, 0.0},
    {"GTG", "Val", 20, {2,3,2}, 0, 0.0},


    // previous talbes use the end codons as different sets. Here, we use them as synonimus codons
    {"TGA", "END", 21, {3,2,0}, 0, 0.0},
    {"TAG", "END", 22, {3,0,2}, 0, 0.0},
    {"TAA", "END", 23, {3,0,0}, 0, 0.0},

};

int main(void) {

    //printf("A");

    char *sequence = NULL;
    size_t seq_capacity = 0;
    size_t seq_length = 0;

    char line[1024];
    while (fgets(line, sizeof(line), stdin)) {
        // skip the header
        if (line[0] == '>') {
            continue;
        }
       
        // remove trailing whitespace/newlines
        char *p = line;
        while (*p) {
            if (*p == '\n' || *p == '\r') {
                *p = '\0';
                break;
            }
            p++;
        }
        // append to  buffer
        size_t line_len = strlen(line);
        if (seq_length + line_len + 1 > seq_capacity) {
            // make it larger so we can keep appending sequences
            seq_capacity = (seq_capacity + line_len + 1024) * 2;
            sequence = realloc(sequence, seq_capacity);
            if (!sequence) {
                fprintf(stderr, "Memory allocation error\n");
                return 1;
            }
        }
        // adjust length
        memcpy(sequence + seq_length, line, line_len);
        seq_length += line_len;
        sequence[seq_length] = '\0';
    }

    // count codons by 3

    long total_codons = 0;
    for (size_t i = 0; i + 2 < seq_length; i += 3) {
        int c0 = base_to_code(sequence[i]);
        int c1 = base_to_code(sequence[i+1]);
        int c2 = base_to_code(sequence[i+2]);
        if (c0 < 0 || c1 < 0 || c2 < 0) {
            //  invalid base => skip and not counting it
            continue;
        }

        // al to upper case
        char codon_str[4];
        codon_str[0] = toupper((unsigned char)sequence[i]);
        codon_str[1] = toupper((unsigned char)sequence[i+1]);
        codon_str[2] = toupper((unsigned char)sequence[i+2]);
        codon_str[3] = '\0';

        // find it in codon_table
        for (int c=0; c<NUM_CODONS; c++) {
            if (strcmp(codon_str, codon_table[c].codon) == 0) {
                codon_table[c].count++;
                total_codons++;
                break;
            }
        }
    }

    // compute frequency for each codon
    for (int i=0; i<NUM_CODONS; i++) {
        if (codon_table[i].count > 0) {
            // freq is calculated by 1000
            codon_table[i].frequency = (codon_table[i].count / (double)total_codons) * 1000.0;
        } else {
            codon_table[i].frequency = 0.0;
        }
    }

    // find the codons with the same aa to get the min, max and the avg frequency
    double min_freq[NUM_CODONS];
    double max_freq[NUM_CODONS];
    double avg_freq[NUM_CODONS];
    int    synonym_count[NUM_CODONS];

    for (int i=0; i<NUM_CODONS; i++) {
        min_freq[i] = 9999999.0;
        max_freq[i] = -9999999.0;
        avg_freq[i] = 0.0;
        synonym_count[i] = 0;
    }

    for (int i=0; i<NUM_CODONS; i++) {
        int aa = codon_table[i].aa_number;
        for (int j=0; j<NUM_CODONS; j++) {
            if (codon_table[j].aa_number == aa) {
                double f = codon_table[j].frequency;
                if (f < min_freq[i]) min_freq[i] = f;
                if (f > max_freq[i]) max_freq[i] = f;
                avg_freq[i] += f;
                synonym_count[i]++;
            }
        }
    }
    //  average
    for (int i=0; i<NUM_CODONS; i++) {
        if (synonym_count[i] > 0) {
            avg_freq[i] /= synonym_count[i];
        } else {
            min_freq[i] = 0.0;
            max_freq[i] = 0.0;
            avg_freq[i] = 0.0;
        }
    }

    // final table
    for (int i=0; i<NUM_CODONS; i++) {

        printf("%d  %d  %d  %.2f  %.2f  %.2f  %.2f  %d  %s  %s\n",
               codon_table[i].base_code[0],
               codon_table[i].base_code[1],
               codon_table[i].base_code[2],
               codon_table[i].frequency,
               max_freq[i],
               min_freq[i],
               avg_freq[i],
               codon_table[i].aa_number,
               codon_table[i].aa_name,
               codon_table[i].codon);
    }

    free(sequence);
    return 0;
}
