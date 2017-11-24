#define FILEOP_C
#include "fileop.h"
#include <string.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_rng.h>

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

void build_filename(char* name, char* directory, char* extension, double num)
{
    char buffer[200];
    snprintf(buffer, 200, "%.2lf", num);
    strcpy(name, directory);
    strcat(name, buffer);
    strcat(name, extension);
}

char* generate_code(int n, gsl_rng* r)
{
    char* code = malloc((n+1)*sizeof(char));
    for(int i=0; i<n; i++)
    {
        int c = 65+26.*gsl_rng_uniform(r);
        code[i] = c;
    }
    code[n] = '\0';

    return code;
}

char* alloc_coded_filename(char* suffix, char* code)
{
    int n = strlen(code);
    int m = strlen(suffix);
    char* name = malloc((n+m+6)*sizeof(char));

    for(int i=0; i<n; i++)
        name[i] = code[i];
    name[n] = '_';
    for(int i=0; i<m; i++)
        name[i+n+1] = suffix[i];
    name[n+m+1] = '.';
    name[n+m+2] = 't';
    name[n+m+3] = 'x';
    name[n+m+4] = 't';
    name[n+m+5] = '\0';

    return name;
}

char* alloc_folder_filename(char* suffix, char* folder)
{
    int n = strlen(folder);
    int m = strlen(suffix);
    char* name;

    // check if / was included in the last folder name
    int check = 0;
    if(folder[n-1] == '/')
        check = 1;

    // allocate memory 
    if(check)
        name = malloc((n+m+1)*sizeof(char));
    else
        name = malloc((n+m+2)*sizeof(char));

    // put folder name
    for(int i=0; i<n; i++)
        name[i] = folder[i];
    
    // if / was included, put suffix
    if(check)
    {
        for(int i=0; i<m; i++)
            name[i+n] = suffix[i];
        name[n+m] = '\0';
    }
    // if / was not included, write it and then write suffix
    else
    {
        name[n] = '/';
        for(int i=0; i<m; i++)
            name[i+n+1] = suffix[i];
        name[n+m+1] = '\0';
    }

    return name;
}



// computes average of columns and writes it in file output
void avgfile(int row, int col, char *input, char *output)
{
    int nr=0;

    FILE *infile;
    FILE *outfile;

    //allocate vector which stores the average value of
    //the columns and initialize its elements to zero
    double *v;
    v = calloc(col, sizeof(double));

    //open file and store the sum of all
    //the elements in a column
    infile = fopen(input, "r");
    for(int i=0; i<row; i++)
    {
        for(int j=0; j<col; j++)
        {
            double temp;
            nr += fscanf(infile, "%lf", &temp);
            v[j] += temp;
        }
    }
    fclose(infile);

    //divide by the number of elements
    for(int j=0; j<col; j++)
        v[j] /= (double)row;


    // write averages on file
    outfile = fopen(output, "w");
    for(int j=0; j<col; j++)
        fprintf(outfile, "%lf\n", v[j]);

    fclose(outfile);
    free(v);

}


// outputs histogram from file
// inname: input file name
// outname: output file name
// Nbin: number of bins in the histogram
// min, max: histogram extrema
// row: number of rows in input file
// col: number of columns in input file
void histo(char* inname, char* outname, int Nbin, int min, int max, int row, int col)
{
    int nr=0;
    
    int* histo;    //histogram
    histo = calloc(Nbin, sizeof(int));
    
    FILE *infile;
    FILE *outfile;
    

    // READ FROM FILE AND GENERATE HISTOGRAM
    
    infile = fopen(inname, "r");
    
    for(int i=0; i<row; i++)
    {
        for(int j=0; j<col; j++)
        {
            double val;
            nr += fscanf(infile, "%lf", &val);
            int idx = (val-min) / ((double)(max-min)/(double)Nbin);
            if(idx < Nbin && idx >= 0)
                histo[idx]++;
        }
    }

    fclose(infile);

    
    // WRITE HISTOGRAM ON FILE
    
    outfile = fopen(outname, "w");

    for(int i=0; i<Nbin; i++)
        fprintf(outfile, "%lf %lf\n", ((max-min)*(double)i/(double)Nbin)+min, (double)histo[i]/(row*col));

    fclose(outfile);
    

    free(histo);
}

// computes the eigenvalues of D from the eigenvalues
// of H in the type(1,0) geometry and writes them on file
void DfromH_t10(int row, int col, char *input, char *output)
{
    int nr=0;

    FILE *infile;
    FILE *outfile;

    //allocate vector which stores the H eigenvalues
    double *v;
    v = malloc(col*sizeof(double));
    //allocate vector which stores the D eigenvalues
    double *w;
    w = malloc(col*col*sizeof(double));

    //open input and output file
    infile = fopen(input, "r");
    outfile = fopen(output, "w");

    //read H eigenvalues, compute D eigenvalues, and write them on file
    for(int i=0; i<row; i++)
    {
        //read H eigenvalues
        for(int j=0; j<col; j++)
        {
            double temp;
            nr += fscanf(infile, "%lf", &temp);
            v[j] = temp;
        }


        //compute D eigenvalues

        //off diagonal
        for(int j=0; j<col; j++)
        {
            for(int k=0; k<j; k++)
            {
                int idx1 = j*col + k;
                int idx2 = k*col + j;
                double temp = v[j] + v[k];
                w[idx1] = temp;
                w[idx2] = temp;
            }
        }

        //diagonal
        for(int j=0; j<col; j++)
        {
            int idx = j*col + j;
            w[idx] = 2.*v[j];
        }


        //order D eigenvalues
        gsl_sort(w, 1, col*col);


        //write D eigenvalues on file
        for(int j=0; j<col*col; j++)
            fprintf(outfile, "%lf ", w[j]);
        fprintf(outfile, "\n");
    }

    
    fclose(infile);
    fclose(outfile);
    free(v);
    free(w);
}






