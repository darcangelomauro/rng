#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "statistics.h"
#include "global.h"

#define N 501

int main(int argc, char** argv)
{
    if(argc != 2)
    {
        printf("gimme the code\n");
        exit(EXIT_FAILURE);
    }

    char* data;
    char* simM;
    int n = strlen(argv[1]);
    data = malloc((n+10)*sizeof(char));
    simM = malloc((n+10)*sizeof(char));
    alloc_input_filename(data, simM, argv[1]);

    FILE* fsimM = fopen(simM, "r");
    if(fsimM == NULL)
    {
        printf("Error, failed to open input file\n");
        exit(EXIT_FAILURE);
    }
    FILE* fout = fopen("out.txt", "w");
    if(fout == NULL)
    {
        printf("Error, failed to open output file\n");
        exit(EXIT_FAILURE);
    }

    read_input(fsimM, 8, 2, 6, N, fout);

    fclose(fsimM);
    fclose(fout);
}
