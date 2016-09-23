#include <stdio.h>
#include <stdlib.h>
#include <string.h>

long double mHG_pvalue (int B, int N, int max_size, double min_hgt)
{
	int N_ = (N < max_size) ? N : max_size;
	int B_ = (B < max_size) ? B : max_size;


	// we allocate m on the heap rather than on the stack, to
	// prevent any memory issue in case B and/or N are large
	long double **m;
	m = (long double **)malloc(sizeof(long double *) * (B_ + 1));
	int i;
	for (i = 0; i < (B_ + 1); i++)
	{
		m[i] = (long double *)malloc(sizeof(long double) * (N_ + 1));
		memset(m[i], 0, sizeof(long double) * (N_ + 1));
	}
	m[0][0] = 1;
		
	long double base_hg = 1;
	
	int n;
	for (n = 1; n <= N_; n++)
	{
		int min_nB;
		if (B >= n)
		{
			min_nB = n;
			base_hg = base_hg * (B - n + 1) / (N - n + 1);
		}
		else
		{
			min_nB = B;
			base_hg = base_hg * n / (n - B);
		}
		
		long double tail_hg = base_hg;
		long double curr_hg = base_hg;

		// first loop - sum up the tail, until the sum is bigger than min_hgt
		int b;
		for (b = min_nB; (tail_hg <= min_hgt) && (b > 0); b--)
		{
			m[b][n] = 0;

			curr_hg = curr_hg * (b * (N - B - n + b)) / ((n - b + 1) * (B - b + 1));
			tail_hg += curr_hg;
		}

		// second loop, starts when b is the maximal for which
		// HGT(N,B,n,b) > min_hgt
		for (; b > 0; b--)
		{
			m[b][n] = 0;

			// calculate current cell value by two optional
			// cells from which it can be reached

			// 1. last element in vector is 0
			if (m[b][n - 1] <= 1)
				m[b][n] += m[b][n - 1] * (N - B - n + b + 1) / (N - n + 1);

			// 2. last element in vector is 1
			if (m[b - 1][n - 1] <= 1)
				m[b][n] += m[b - 1][n - 1] * (B - b + 1) / (N - n + 1);
		}

		m[b][n] = 0;
		m[0][n] += m[0][n - 1] * (N - B - n + 1) / (N - n + 1);
	}

	long double r = 0, temp = 0, y = 0, e = 0;
	for (i = 0; i < B_ + 1; i++) {
		temp = r;
		y = m[i][N_] + e;
		r = temp + y;
		e = (temp - r) + y;
		//r += m[i][N_];
	}

	
	for (i = 0; i < (B_ + 1); i++)
	{
		free(m[i]);
		m[i] = NULL;
	}
	free(m);
	m = NULL;
	
	
	return (1.0 - r);
}

int main(int argc, char *argv[]) {
	if(argc != 4) {
		fprintf(stderr, "Syntax: ./mHG_pvalue mHG_score N B\n");
		exit(-1);
	}
	
	double mHG_score = atof(argv[1]);
	int N = atoi(argv[2]);
	int B = atoi(argv[3]);

	printf("pval = %e\n", mHG_score);
	
	long double pval = mHG_pvalue (B, N, N, mHG_score);
	printf("%.6Le\n", pval);
	
	
	return 0;
}
