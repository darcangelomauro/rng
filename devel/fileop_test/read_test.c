#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fileop.h"


int main(int argc, char** argv)
{

    if(argc < 2)
    {
        printf("gimme code\n");
        exit(EXIT_FAILURE);
    }

    char* file = alloc_coded_filename("varG_args", argv[1]);
    FILE* f = fopen(file, "r");
    char buffer[200];

    int j = 0;
    while(fscanf(f, "%s", buffer) != EOF)
    {
        if(j > 40) break;
        printf("code: %s, length: %d\n", buffer, (int)strlen(buffer));
        j++;
    }

}
