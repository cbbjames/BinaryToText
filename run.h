// Paul Draghicescu
// pdraghicescu@pnri.org
//

#ifndef RUN_H 
#define RUN_H

#define TUP_BUF_SIZE 10000

typedef struct option {
    int n, k;
    double b, s_std;
    int s_option;
    char in_file1[50];
    char in_file2[50];
} option_t;

int get_header(FILE* fbin,  int *d, int *n_vars, unsigned long *n_tuples, double *average);
int run(option_t *opt);
double calculate_std_bin(FILE* fbin, double average, int d, unsigned long n_tuples);

#endif

