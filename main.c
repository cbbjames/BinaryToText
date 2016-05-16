// main.c
//
// Paul Draghicescu
// pdraghicescu@pnri.org
// 4/6/2015
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "run.h"
FILE *input1, *input2;
int d, n_vars;
unsigned long n_tups;
double avg;
/*
void read_input2() {
	char* mbuf;
	long int len = lseek(fd, 0, SEEK_END);
	mbuf = (char*)mmap(NULL, len, PROT_READ, MAP_PRIVATE, input2, 0);
	analyse(mbuf, len);
}
*/
int parse_args(int argc, char* argv[], option_t* opt) { 
    int c;
    while (--argc > 0 && (*++argv)[0] == '-')
        while ((c = *++argv[0]))
            switch (c) {
                case 'n':
                    {
                        if (--argc > 0) {
                            opt->n = (atoi(argv[1]));
                            argv++;
                            *argv+= strlen(*argv)-1;
                        }
                        break;
                    }
                case 'k':
                    {
                        if (--argc > 0) {
                            opt->k = (atoi(argv[1]));
                            argv++;
                            *argv+= strlen(*argv)-1;
                        }
                        break;
                    }
                case 'b':
                    {
                        if (--argc > 0) {
                            opt->b = (atof(argv[1]));
                            argv++;
                            *argv+= strlen(*argv)-1;
                        }
                        break;
                    }
                case 's':
                    {
                        if (--argc > 0) {
                            opt->s_option = 1;
                            opt->s_std = (atof(argv[1]));
                            argv++;
                            *argv+= strlen(*argv)-1;
                        }
                        break;
                    }
                default:
                    {
                        fprintf(stderr, "Illegal option\n");
                        argc = 0;
                        break;
                    }
            }
    if (argc != 2) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   run [options] <input1> <input2>\n\n");
        fprintf(stderr, "Options:                                                 \n");
        fprintf(stderr, "         -n INT    number of sorted tuples for out1.txt and out2.txt [%d]\n", opt->n);
        fprintf(stderr, "         -k INT    number of sorted tuples for out3.txt [%d]\n", opt->k);
        fprintf(stderr, "         -b FLOAT  bin width for out4.txt [%.2f]\n", opt->b);
        fprintf(stderr, "         -s FLOAT  number of standard deviations for out5.txt[%.2f]\n", opt->s_std);
        return -1;
    } 
    strcpy(opt->in_file1, argv[0]);
    strcpy(opt->in_file2, argv[1]);
    return 0;
}
                    
int init_options(option_t *opt) {
    opt->n = 0;
    opt->k = 0;
    opt->b = 0;
    opt->s_std     = 0.0;
    opt->s_option  = 0;
    return 0;
}

int main(int argc, char * argv[]) {
    int err;
    time_t ltime;
    option_t opt;//n,k,b,s

    err = init_options(&opt);
    if (err < 0) {
        fprintf(stderr, "Error initializing arguments.\n");
        return 1;
    }

    err = parse_args(argc, argv, &opt);
    if (err < 0) {
        return 1;
    }

    //TODO: remove this, implement these options
    if ((opt.n != 0) || (opt.k != 0) || (opt.b != 0)) {
        fprintf(stderr, "WARNING: -n, -k, -b not yet implemented!\n");
        //return 1;
    }
    ltime = time(NULL); 
    fprintf(stderr, "%s",asctime(localtime(&ltime)));

    err = run(&opt);
    if (err < 0) {
        fprintf(stderr, "Error running.\n");
        return 1;
    }
    fprintf(stderr, "Finished Succesfully.\n");
    return 0;
}

