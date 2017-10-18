#define DISKMAP_C

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>

void *load_mmap(char const *filename, int *fd, size_t size, int make_room)
{
    *fd = open(filename, make_room ? O_RDWR | O_CREAT | O_TRUNC : O_RDWR, (mode_t)0600);
    if(*fd==-1)
    {
        printf("Error opening file\n");
        return NULL;
    }

    if(make_room)
    {
        int result=lseek(*fd, size-1, SEEK_SET);
        if(*fd==-1)
        {
            printf("Error stretching file\n");
            close(*fd);
            return NULL;
        }

        result = write(*fd, "", 1);
        if(result!=1)
        {
            printf("Error writing file\n");
            close(*fd);
            return NULL;
        }
    }

    void *map=mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, *fd, 0);
    if(map==MAP_FAILED)
    {
        printf("Error mapping file\n");
        return NULL;
    }

    return map;
}

int release_mmap(void *map, size_t size, int fd)
{
    if(munmap(map, size) == -1)
    {
        printf("Error un-mapping file\n");
        return -1;
    }

    close(fd);

    return 0;
}

// creates readable file from mmapped array.
// can also organize data in rows with multiple columns.
// just be sure that rows*cols == array lenght, otherwise
// bad things are gonna happen
void readable(double map[], char const *filename, int rows, int cols)
{
    FILE* frble;
    frble = fopen(filename, "w");

    for(int i=0; i<rows; i++)
    {
        for(int j=0; j<cols; j++)
        {
            fprintf(frble, "%.17lf ", map[i*cols + j]);
        }

        fprintf(frble, "\n");
    }

    fclose(frble);
}
                    
