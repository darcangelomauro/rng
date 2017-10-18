#ifndef DATA_H
#define DATA_H

extern void print_data(FILE* fdata);
extern void print_time(FILE* fdata, char* s);
extern void print_thermalization(FILE* fobsS);
extern void print_thermalization_plus(FILE* fobsStr);
extern void print_simulation(FILE* fobsS, FILE* fobsHL);
#endif
