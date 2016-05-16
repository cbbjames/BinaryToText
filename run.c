// Paul Draghicescu
// pdraghicescu@pnri.org
//

#include <math.h>
#include <stdio.h>
#include <float.h>
#include "run.h"
#include "string.h"
#include<stdlib.h>

#ifdef DELTA_PREC
#define PRINT_PRC_D "%.10Lf%c"
#else
#define PRINT_PRC_D "%.10lf%c"
#endif

// write all tuples above cutoff
unsigned long count_tuples_bin_cutoff(const option_t* opt, FILE* fbin, int d, double cutoff, long unsigned n_tuples, double average, double std) {
    double next;
    unsigned long c = 0;
    int skip = d+sizeof(double)/sizeof(int);
    int i,j;
    int *ip;
    int bufsize = (d*TUP_BUF_SIZE)+((sizeof(double)/sizeof(int))*TUP_BUF_SIZE);
    int buffer[bufsize];
    unsigned int n;
    FILE* fout5 = fopen("out5.txt","w");
    if (!fout5) {
        fprintf(stderr, "Error: couldn't open out5.txt\n");
        return -1;
    }
    while(n_tuples > 0)  {
        n = fread(buffer, sizeof(int), bufsize,  fbin);
        if (!n) break;
        for (j = 0; (j < n) && (--n_tuples > 0); j+=skip) {
            memcpy(&next,&buffer[j+d],sizeof(double));
            // if next score is in the tail
            if ((!(d%2) && (next >= cutoff)) || ((d%2) && (next <= cutoff))) {
                //TODO: inefficient, need different types for different d
                for (i = 0, ip = &buffer[j]; i < d; i++, ip++) {
                    fprintf(fout5,"%d\t",*ip);
                }
                fprintf(fout5, PRINT_PRC_D,next,'\t');
                fprintf(fout5, PRINT_PRC_D,fabs(next-average)/std, '\n');
                c++;
            }
        }
    }
    fclose(fout5);
    return c;
}
// calculate standard deviation, min and max delta
double calculate_std_bin(FILE* fbin, double average, int d, unsigned long n_tuples) {
	double next;
	unsigned long c = 0;
	double sum = 0.0;
	double result;
	double min = DBL_MAX;
	double max = -DBL_MAX;
	int skip = d + sizeof(double) / sizeof(int);
	int j;
	// TODO: hard code types for different values of d
	int bufsize = (d*TUP_BUF_SIZE) + ((sizeof(double) / sizeof(int))*TUP_BUF_SIZE);
	int *buffer;
	buffer = (int*)malloc(bufsize * sizeof(int));
	unsigned int n;
	while (n_tuples > 0) {
		n = fread(buffer, sizeof(int), bufsize, fbin);
		if (!n) break;
		for (j = 0; (j < n) && (--n_tuples > 0); j += skip) {
			memcpy(&next, &buffer[j + d], sizeof(double));
			if (next < min)
				min = next;
			if (next > max)
				max = next;
			sum += pow(next - average, 2);
			c++;
		}
	}
	fprintf(stderr, "Min Delta:\t\t");
	fprintf(stderr, PRINT_PRC_D, min, '\n');
	fprintf(stderr, "Max Delta:\t\t");
	fprintf(stderr, PRINT_PRC_D, max, '\n');
	result = sum / c;
	free(buffer);
	return sqrtl(result);
}
// read and store header information
int get_header(FILE* fbin, int *d, int *n_vars, unsigned long *n_tuples, double *average) {
    int err;
    err = fread(d, 1, sizeof(int), fbin);
    if (!err) {
        fprintf(stderr, "Error: couldn't read dummy 0 from file\n");
        return -1;
    }
    err = fread(d, 1, sizeof(int), fbin);
    if (!err) {
        fprintf(stderr, "Error: couldn't read dimension from file\n");
        return -1;
    }
    err = fread(n_vars, 1, sizeof(int), fbin);
    if (!err) {
        fprintf(stderr, "Error: couldn't read n_vars from file\n");
        return -1;
    }
    err = fread(n_tuples, 1, sizeof(unsigned long), fbin);
    if (!err) {
        fprintf(stderr, "Error: couldn't read n_tuples from file\n");
        return -1;
    }
    err = fread(average, 1, sizeof(double), fbin);
    if (!err) {
        fprintf(stderr, "Error: couldn't read average from file\n");
        return -1;
    }
    return 0;
}
int cmp_1(const void*a,const void*b) {
	return (*(double**)b - *(double**)a);
}//big->small
void analyse(char *buf, int len,int* numbers)
{
	int i;
	numbers[i = 0] = 0;
	for (char *p = buf; *p && p - buf<len; p++)
		if (*p == ' ')
			numbers[++i] = 0;
		else
			numbers[i] = numbers[i] * 10 + *p - '0';
}
char *fastStrcat(char *pszDest, const char* pszSrc)
{
	while (*pszDest)
	      pszDest++;
	while ((*pszDest++ = *pszSrc++));
		return --pszDest;
}
void wbuf(char* buf, int *tup, int d) {
	char buff[64];
	for (int i = 0; i < d; i++) {
		sprintf(buff, "%d", tup[i]);
		buf = fastStrcat(buf, buff);
	}
	buf = fastStrcat(buf, "\n");
}
void out1(FILE* fout1,int n, double* *keyvalue, double *tscore,int d, int **readin) {
	for (int i = 0; i < n; i++) {
		char*buf;
		buf = (char*)malloc(3 * (d + 5) * sizeof(char));
		buf[0] = '\0';
		int index = keyvalue[i] - &tscore[0];
		char t[2];
		for (int j = 0; j < d - 1; j++) {
			sprintf(t, "%d\t", readin[index][j]);
			buf = fastStrcat(buf, t);
		}
		sprintf(t, "%d\n", readin[index][d - 1]);
		buf = fastStrcat(buf, t);
		fwrite(buf, sizeof(char), strlen(buf), fout1);
		free(buf);
	}
}
void out2(FILE* fout2, int n, double* *keyvalue, double *tscore, int d, int **readin,int n_tuples) {
	for (int i = n_tuples - 1; i >= n_tuples - n; --i) {
		char*buf;
		buf = (char*)malloc(3 * (d + 5) * sizeof(char));
		buf[0] = '\0';
		int index = keyvalue[i] - &tscore[0];
		char t[2];
		for (int j = 0; j < d - 1; j++) {
			sprintf(t, "%d\t", readin[index][j]);
			buf = fastStrcat(buf, t);
		}
		sprintf(t, "%d\n", readin[index][d - 1]);
		buf = fastStrcat(buf, t);
		fwrite(buf, sizeof(char), strlen(buf), fout2);
		free(buf);
	}
}
void out3(FILE* fout3, int n, double* *keyvalue, double *tscore, 
	int d, int **readin, int n_tuples,option_t *opt) {
	int total_cnt = 0;
	double avg3 = 0;
	int k = opt->k, k2 = opt->k;
	char *ans;
	ans = (char*)malloc(sizeof(char)*(2 * k + 1)*d);
	ans[0] = 0;
	int index;
	//max side
	for (int i = 0; i < n_tuples; i++) {
		if (!k) break;
		index = keyvalue[i] - &tscore[0];
		for (int j = 0; j < d; j++)
			if (readin[index][j] == k) { k--; total_cnt++; break; }//find
		wbuf(ans, readin[index], d);
		avg3 += *keyvalue[i];
	}
	//min side
	for (int i = n_tuples - 1; i >= 0; --i) {
		if (!k2) break;
		index = keyvalue[i] - &tscore[0];
		for (int j = 0; j < d; j++)
			if (readin[index][j] == k) { k--; total_cnt++; break; }//find
		wbuf(ans, readin[index], d);
		avg3 += *keyvalue[i];
	}
	avg3 /= total_cnt;
	char first_line[64];
	sprintf(first_line, "%d\t%f", k, avg3);
	strcat(first_line, ans);
	fwrite(first_line, sizeof(char), strlen(first_line), fout3);
}
void out4(FILE* fout4, int n, double* *keyvalue, double *tscore,
	int d, int **readin, int n_tuples, option_t *opt) {

	double max = *keyvalue[0], min = *keyvalue[n_tuples - 1];
	double b = opt->b;
	int line_num = floor((max - min) / b) + 1;
	double **min_max;//[min,max)
	min_max = (double**)malloc(line_num * sizeof(double*));
	for (int i = 0; i < line_num; i++) {//0=min 1=max
		min_max[i] = (double*)malloc(2 * sizeof(double));
		min_max[i][0] = min + i*b;
		min_max[i][1] = min + (i +1)* b;
	}

	int each_cnt[line_num];
	int total_cnt = 0;
	for (int i = 0; i < line_num; i++) {
		int cnt = 0;
		while (*keyvalue[total_cnt] < min_max[i][1]) {
			cnt++;
			total_cnt++;
		}
		each_cnt[i] = cnt;
	}

	char first_line[64];
	sprintf(first_line, "%f\t%f\n", min, max);
	fwrite(first_line, sizeof(char), strlen(first_line), fout4);
	char ans[64 * line_num];
	ans[0] = 0;
	for (int i = 0; i < line_num; i++) {
		sprintf(first_line, "%d\n", each_cnt[i]);
		strcat(ans, first_line);
	}
	fwrite(ans, sizeof(char), strlen(ans), fout4);
}
void out(FILE* fbin2,double average, int d, unsigned long n_tuples,int n, option_t *opt) {

	double *tscore;
	double* *keyvalue;//score index->index
	tscore = (double*)malloc(n_tuples * sizeof(double));
	keyvalue = (double**)malloc(n_tuples * sizeof(double*));
	int **readin;
	readin = (int**)malloc(n_tuples * sizeof(int*));
	for (int i = 0; i < n_tuples; i++)
		readin[i] = (int*)malloc((d + 5) * sizeof(int));
	for (int i = 0; i < n_tuples; i++) {
		char *buf;
		buf = (char*)malloc((d + 5) * sizeof(char));
		fgets(buf, d + 5, fbin2);
		analyse(buf, sizeof(buf), readin[i]);
		fscanf(fbin2, "%f", &tscore[i]);
		keyvalue[i] = &tscore[i];
		free(buf);
	}
	qsort((void*)keyvalue, n_tuples, sizeof(double*), cmp_1);
	//------------------------------------------------------------------------
	FILE* fout1 = fopen("out1.txt", "w");
	FILE* fout2 = fopen("out2.txt", "w");
	if (n == 0) { fclose(fout1); fclose(fout2); }
	//1---------------------------------------------------------------------
	out1(fout1, n,keyvalue,tscore,d,readin);
	fclose(fout1);
	//2--------------------------------------------------------------------------
	out2(fout2, n, keyvalue, tscore, d, readin,n_tuples);
	fclose(fout2);
	//3--------------------------------------------------------------------------
	FILE* fout3 = fopen("out3.txt", "w");
	out3(fout2, n, keyvalue, tscore, d, readin, n_tuples,opt);
	fclose(fout3);
	//4------------------------------------------------------------------------------
	FILE* fout4 = fopen("out4.txt", "w");
	out4(fout2, n, keyvalue, tscore, d, readin, n_tuples, opt);
	fclose(fout4);
	//
	free(tscore); free(readin); free(keyvalue);
}

// main function 
int run(option_t *opt) {
    FILE* fbin1;
    FILE* fbin2;
    double average = 0.0;
    double std = 0.0;
    double cutoff;
    int d, n_vars, sign=1;
    unsigned long n_tuples, c;

    // open output file, get dimension and other information from header
    fbin1 = fopen(opt->in_file1, "rb");
    if (!fbin1) {
        fprintf(stderr, "Error: couldn't open header file \"%s\".\n",opt->in_file1);
        return -1;
    }
    fbin2 = fopen(opt->in_file2, "rb");
    if (!fbin2) {   
        fprintf(stderr, "Error: couldn't open binary file 2\"%s\".\n",opt->in_file1);
        return -1;
    }
    fseek(fbin1, 0L, SEEK_END);
    if (!ftell(fbin1)) {
        fprintf(stderr, "Error: empty binary file 1  \n");
        return -1;
    }
    fseek(fbin2, 0L, SEEK_END);
    if (!ftell(fbin2)) {
        fprintf(stderr, "Error: empty binary file 2\n");
        return -1;
    }
    rewind(fbin1);
    rewind(fbin2);
    get_header(fbin1, &d, &n_vars, &n_tuples, &average);
    if (fclose(fbin1) != 0) {
        fprintf(stderr, "Error closing binary file 2.\n");
        return -1;
    }
    sign = pow(-1,d);
    fprintf(stderr, "Dimension:\t\t%d\n",d);
    fprintf(stderr, "Number of Variables:\t%d\n",n_vars);
    fprintf(stderr, "Number of Tuples:\t%ld\n",n_tuples);
    fprintf(stderr, "Average:\t\t%.10f\n",average);

	out(fbin2,average,d,n_tuples,opt->n,opt);

	//sample 5
    // if -s is passed, pass through the file once to calculate std, and again to create out5.txt
    if (opt->s_option) {
        std = calculate_std_bin(fbin2, average, d, n_tuples+1);
        rewind(fbin2);
        cutoff = average + sign*std*opt->s_std;
        fprintf(stderr, "St. Dev.:\t\t%.10f\n",std);
        fprintf(stderr, "Cutoff:\t\t\t%.10f\n",cutoff);

        // count and print tuples above cutoff
        c = count_tuples_bin_cutoff(opt, fbin2, d, cutoff, n_tuples+1, average,std);
        fprintf(stderr, "Tuples Above Cutoff:\t%ld\t(%.2f%%)\n", c, 100.0*c/n_tuples);
    }

    if (fclose(fbin2) != 0) {
        fprintf(stderr, "Error closing binary file 2.\n");
        return -1;
    }
    return 0;
}

