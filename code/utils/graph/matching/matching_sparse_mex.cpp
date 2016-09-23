#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "mex.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
typedef int mwSize;
#endif // MX_API_VER


/**
 * n the number of nodes
 * m the number of nodes
 * nedges the number of edges
 * vv1 is the source for each of the nedges 
 * vv2 is the target for each of the nedges
 * weight is the weight of each of the nedges
 * out1 is a vector of length at most min(n,m),
 * out2 is a vector of length at most min(n,m),
 * noutedges is the number of out edges
 */
double match(int n, int m, int nedges, double *vv1, double *vv2, double *weight, double *out1, double *out2,
             int *noutedges)
{
	double ret, al;
	double *l1, *l2, *w;
	int *match1, *match2, *v1, *v2;
	int i, j, k, p, q, r, t1, t2;
	int *s, *t, *deg, *offset, *list;
	bool *tight;

	l1 = new double[n];
	l2 = new double[n+m];
	v1 = new int[nedges];
	v2 = new int[nedges];
	s = new int[n];
	t = new int[n+m];
	match1 = new int[n];
	match2 = new int[n+m];
	offset = new int[n];
	deg = new int[n];
	list = new int[nedges + n];
	w = new double[nedges + n];
	tight = new bool[nedges + n];

	for (i = 0; i < nedges; i++) {
		v1[i] = (int)(vv1[i] + .5);
		v2[i] = (int)(vv2[i] + .5);
	}
	for (i = 0; i < n; i++) {
		offset[i] = 0;
		deg[i] = 1;
	}
	for (i = 0; i < nedges; i++) deg[v1[i]]++;
	for (i = 1; i < n; i++) offset[i] = offset[i-1] + deg[i-1];
	for (i = 0; i < n; i++) deg[i] = 0;
	for (i = 0; i < nedges; i++) {
		list[offset[v1[i]] + deg[v1[i]]] = v2[i];
		w[offset[v1[i]] + deg[v1[i]]] = weight[i];
		deg[(int)v1[i]]++;
	}
	for (i = 0; i < n; i++) {
		list[offset[i] + deg[i]] = m + i;
		w[offset[i] + deg[i]] = 0;
		deg[i]++;
	}
	for (i = 0; i < n; i++) {
		l1[i] = 0;
		for (j = 0; j < deg[i]; j++) {
			if (w[offset[i]+j] > l1[i]) l1[i] = w[offset[i] + j];
		}
	}
	for (i = 0; i < n; i++) {
		match1[i] = -1;
	}
	for (i = 0; i < n + m; i++) {
		l2[i] = 0;
		match2[i] = -1;
	}
	for (i = 0; i < n; i++) {
		for(j = 0; j < n + m; j++) t[j] = -1;
		s[p = q = 0] = i;
		for(; p <= q; p++) {
			if (match1[i] >= 0) break;
			k = s[p];
			for (r = 0; r < deg[k]; r++) {
				if (match1[i] >= 0) break;
				j = list[offset[k] + r];
				if (w[offset[k] + r] < l1[k] + l2[j] - 1e-8) continue;
				if (t[j] < 0) {
					s[++q] = match2[j];
					t[j] = k;
					if (match2[j] < 0) {
						for(; j>=0 ;) {
							k = match2[j] = t[j];
							p = match1[k];
							match1[k] = j;
							j = p;
						}
					}
				}
			}
		}
		if (match1[i] < 0) {
			al = 1e20;
			for (j = 0; j < p; j++) {
				t1 = s[j];
				for (k = 0; k < deg[t1]; k++) {
					t2 = list[offset[t1] + k];
					if (t[t2] < 0 && l1[t1] + l2[t2] - w[offset[t1] + k] < al) {
						al = l1[t1] + l2[t2] - w[offset[t1] + k];
					}
				}
			}
			for (j = 0; j < p; j++) l1[s[j]] -= al;
			for (j = 0; j < n + m; j++) if (t[j] >= 0) l2[j] += al;
			i--;
			continue;
		}
	}
	ret = 0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < deg[i]; j++) {
			if (list[offset[i] + j] == match1[i]) {
				ret += w[offset[i] + j];
			}
		}
	}
    *noutedges = 0;
    for (i = 0; i < n; i++) {
        if (match1[i] < m) (*noutedges)++;
    }
    *noutedges = 0;
    for (i = 0; i < n; i++) {
        if (match1[i] < m) {
            out1[*noutedges] = i;
            out2[*noutedges] = match1[i];
            (*noutedges)++;
        }
    }
	return ret;
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    int n, m, nedges;
    double *vv1, *vv2, *w, *out1, *out2;
    int noutedges;
    
    int curarg = 0;
    
    n = (int)mxGetScalar(prhs[curarg++]);
    m = (int)mxGetScalar(prhs[curarg++]);
    nedges = (int)mxGetScalar(prhs[curarg++]);
    
    vv1 = mxGetPr(prhs[curarg++]);
    vv2 = mxGetPr(prhs[curarg++]);
    w = mxGetPr(prhs[curarg++]);
    
    int minnm = n;
    if (m < minnm) minnm = m;
    
    plhs[0] = mxCreateDoubleMatrix(minnm,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(minnm,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
    
    double val = match(n,m,nedges,vv1,vv2,w,mxGetPr(plhs[0]),mxGetPr(plhs[1]),&noutedges);
    
    *(mxGetPr(plhs[2])) = (double)noutedges;
    *(mxGetPr(plhs[3])) = val;
}
    

/*
double g[600][400];

int main() {
	int i, j, t=0;
	double k;

	double vv1[10000];
	double vv2[10000];
	double weight[10000];
	double *out1,*out2;

	freopen("data.txt","r",stdin);
	for(i=0;i<600;i++){
		for(j=0;j<400;j++){
			scanf("%lf",&k);
			if (k > 1e-8){
				vv1[t] = i;
				vv2[t] = j;
				weight[t] = k;
				t++;
			}
		}
	}

	printf("%lf\n",match(600,400,t,vv1,vv2,weight,out1,out2));
	return 0;
}
*/

