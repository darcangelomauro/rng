#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "statistics.h"
#include "update.h"
#include "global.h"
#include "macros.h"

int main(int argc, char** argv)
{
    if(argc != 2)
    {
        printf("gimme the code\n");
        exit(EXIT_FAILURE);
    }

    // avoid overwriting previous analysis
    char* checkM = alloc_coded_filename("evalM", argv[1]);
    char* checkD = alloc_coded_filename("evalD", argv[1]);
    FILE* fcheckM = fopen(checkM, "r");
    FILE* fcheckD = fopen(checkD, "r");
    if(fcheckM != NULL || fcheckD != NULL)
    {
        printf("The analysis has already been done\n");
        exit(EXIT_SUCCESS);
    }
    free(checkM);
    free(checkD);
    
    int control[3] = {1,1,0};

    analysis(argv[1], control, P_gamma);

}
