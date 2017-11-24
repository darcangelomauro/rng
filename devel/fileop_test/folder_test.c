#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include "fileop.h"


int main(int argc, char** argv)
{

    if(argc < 3)
    {
        printf("gimme folder\n");
        exit(EXIT_FAILURE);
    }

    for(int i=0; i<5; i++)
    {
        char* code = alloc_folder_filename(argv[1], argv[2]);
        printf("code: %s\n", code);
        char* file = alloc_coded_filename("data", code);
        printf("file: %s\n", file);
        free(code);
        free(file);
    }

}
