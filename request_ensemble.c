#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <curl/curl.h> // guarda que acá en el gcc tengo que pasar manualmente donde esta (path de miniconda)

#define API_URL_FORMAT "https://rest.ensembl.org/sequence/id/%s?object_type=transcript;type=cds;content-type=text/x-fasta"  // dynamic ID URL

// Structure to handle HTTP response
struct MemoryStruct {
    char *memory;
    size_t size;
};

// Callback function for libcurl to write the received data
static size_t WriteMemoryCallback(void *contents, size_t size, size_t nmemb, void *userp) {
    size_t realsize = size * nmemb;
    struct MemoryStruct *mem = (struct MemoryStruct *)userp;

    char *ptr = realloc(mem->memory, mem->size + realsize + 1);
    if(ptr == NULL) {
        // Out of memory
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

    // libcurl
    curl_handle = curl_easy_init();

    if(curl_handle) {
        fprintf(stderr,"Fetching sequence...\n");
        curl_easy_setopt(curl_handle, CURLOPT_URL, url);
        fprintf(stderr,".\n");
        curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, WriteMemoryCallback);
        fprintf(stderr,"..\n");
        curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, (void *)&chunk);
        fprintf(stderr,"...\n");
        curl_easy_setopt(curl_handle, CURLOPT_USERAGENT, "libcurl-agent/1.0");
        fprintf(stderr,"....\n");
       
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

int main(int argc, char *argv[]) {
    if (argc != 2) {  // Expect a single argument for the protein ID
        fprintf(stderr, "Usage: %s <protein_id>\n", argv[0]);
        return 1;
    }

    char *protein_id = argv[1];
    char *sequence_data = fetchProteinSequence(protein_id);  // Pass the protein ID

    if (sequence_data) {
        printf("Fetched Sequence:\n%s\n", sequence_data);
        free(sequence_data);  // Free the allocated memory after use
    } else {
        fprintf(stderr, "Failed to fetch sequence for protein ID: %s\n", protein_id);
        return 1;
    }

    return 0;
}
