/*
 * Filename: buildMST.cpp
 * Author:   Jacky Chen
 * Date:     8/6/2009
 */

/* 
 *	buildMST.cpp :	Finds the maximum weight spanning tree of the graph composed of nodes representing a binary 
 *					feature in the vocabulary and edges representing the mutual information between the two features. 
 */

#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <queue>
#include <algorithm>
#include <vector>
#include <hash_map>
#include "mat.h"
using namespace std;

#pragma comment(lib, "libmat.lib")
#pragma comment(lib, "libmx.lib")

/*
	Computes the mutual information for binary features i and j.
	Input:
		r - # of entries
		ni, nj - # of images having features i, j respectively
		oo - # of cooccurrences of features i and j
*/
double mi(double r, double ni, double nj, double oo) 
{
	double zz, zo, oz, I;

	zz = (r-ni-nj+oo)/r;	//the probability of when both i and j are 0
	zo = (nj-oo)/r;			//...when i is 0 and j is 1
	oz = (ni-oo)/r;			//...when i is 1 and j is 0
	oo /= r;				//...when both i and j are 1

	I = 0;
	if (zz != 0) I += zz*log(zz/((r-ni)/r*(r-nj)/r));
	if (zo != 0) I += zo*log(zo/((r-ni)/r*nj/r));
	if (oz != 0) I += oz*log(oz/(ni/r*(r-nj)/r));
	if (oo != 0) I += oo*log(oo/(ni/r*nj/r));
	I /= log(2.0);
	return I;
}

struct Cost {
	double cost;
	int index;
};

class CompareCost {
public:
	bool operator()(Cost& c1, Cost& c2)
	{
		if (c1.cost < c2.cost) return true;
		return false;
	}
};

/*
	Prim's algorithm on finding the maximum weight spanning tree; returns the parent of each node.
*/
int *prim(int c, double **I){
	double *cost;
	int *prev, *mark;
	cost = new double[c];
	prev = new int[c];
	mark = new int[c];
	int i;
	for (i=0; i<c; i++) {
		cost[i] = 0;
		prev[i] = 0;
		mark[i] = 0;
	}
	cost[0] = 1;
	/* 
		double cost[c] = { 1 };
		int prev[c] = { 0 };
		int mark[c] = { 0 };
	*/

	priority_queue<Cost, vector<Cost>, CompareCost> pq;
	Cost init = {1, 0};
	pq.push(init);
	for (i=1; i<c; i++) {
		Cost init = { 0, i };
		pq.push(init);
	}

	int *tmp;
	int x = 0;
	while (!pq.empty()) {
		if (x < 100 || x%100 == 0) printf("%d edges computed\n", x);
		x++;
		tmp = find(mark, mark+c, 0);
		if (tmp == mark+c) break;

		Cost c1 = pq.top();
		pq.pop();
		int ind = c1.index;
		mark[ind] = 1;
		for (i=0; i<c; i++) {
			if (mark[i] == 1) continue;
			if (i < ind && cost[i] < I[i][ind]) {
				cost[i] = I[i][ind];
				prev[i] = ind;
				Cost c2 = { cost[i], i };
				pq.push(c2);
			} else if (i > ind && cost[i] < I[ind][i]) {
				cost[i] = I[ind][i];
				prev[i] = ind; 
				Cost c2 = { cost[i], i };
				pq.push(c2);
			}
		}
	}

	return prev;
}

/*
	Input: .mat file from convertImagesToObservations
	Output: .mat file specifying parent of each node in the MST, used by fabmap
*/
int main(int argc, char** argv)
{
	string fileName(argv[1]);
	MATFile *mfp = matOpen(fileName.c_str(), "r");
	mxArray *M = matGetVariable(mfp, "observations");

	double *pr;
	int *ip;
	mwSize i, n, j, r, c, k;
	pr = mxGetPr(M);
	n = mxGetNumberOfElements(M);
	r = mxGetM(M);
	c = mxGetN(M);
	ip = (int *)mxMalloc(n * sizeof(*ip));
	for(i=0; i<n; i++) ip[i] = pr[i];

	int *N;
	N = new int[c];
	for(i=0; i<c; i++) {
		N[i] = 0;
		for (j=0; j<r; j++)	N[i] += ip[i*r+j];
	}
	
	//Computing cooccurrences
	int **Nuv;
	Nuv = new int*[c];
	for (i=0; i<c; i++) *(Nuv+i) = new int[c];
	for (i=0; i<c; i++) {
		if (i != 0 && i%100 == 0) printf("Computing cooccurrences up to column %d...\n", i);
		for (j=i+1; j<c; j++) {
			Nuv[i][j] = 0;
			for (k=0; k<r; k++)	Nuv[i][j] += (ip[i*r+k]&&ip[j*r+k]);
		}
	}

	mxFree(ip);
	
	//Computing mutual informations
	int input;
	double **I;
	double ni, nj, oo;
	I = new double*[c];
	for (i=0; i<c; i++) *(I+i) = new double[c];
	for (i=0; i<c; i++) {
		if (i != 0 && i%1000 == 0) printf("Computing mutual informations on features up to %d\n", i);
		for (j=i+1; j<c; j++) {
			ni = N[i];
			nj = N[j];
			oo = Nuv[i][j];
			I[i][j] = mi(r, ni, nj, oo);
		}
	}

	int *prev;
	prev = prim(c, I);

	delete[] I;

	double *dp;
	dp = (double *)mxMalloc(c * sizeof(*dp));
	for(i=0; i<c; i++) 	dp[i] = prev[i];
	
	int pos = fileName.find("obs");
	string dataset = fileName.substr(pos+3);
	string oFileName = "mst";
	mfp = matOpen(oFileName.append(dataset).c_str(), "w");
	mxArray *PREV = mxCreateDoubleMatrix(c, 1, mxREAL);
	double *data = mxGetPr(PREV);
	for (i=0; i<c; i++) data[i] = dp[i];
	matPutVariable(mfp, "prev", PREV);

	printf("Complete.");

	return 0;
}

