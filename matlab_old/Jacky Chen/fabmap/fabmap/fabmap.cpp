/*
 * Filename: fabmap.cpp
 * Author:   Jacky Chen
 * Date:     11/23/2010
 */

/* 
 *	fabmap.cpp :	 Identifies image pairs that form a loop closure using a probabilistic framework and
 *					 recursive Bayesian inference. 
 */

#include "stdafx.h"
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <queue>
#include <algorithm>
#include <vector>
#include <time.h>
#include "mat.h"
using namespace std;

#pragma comment(lib, "libmat.lib")
#pragma comment(lib, "libmx.lib")

void fill_with_unique_rand(int* a, int N, int M)
{
	for(int i=0; i<N; i++) {
		if( (rand() % (N-i)) < M) {// we have to choose required out of the remaining (available-i)
			--M;
			a[M] = i ;
		}
	}
}

struct Entry {
	int obs;
	int loc;
	int count;
};

class CompareEntry {
public:
	bool operator()(Entry& e1, Entry& e2)
	{
		if (e1.count < e2.count) return true;
		return false;
	}
};

int main(int argc, char** argv)
{
	srand(time(NULL));

	string obsFile(argv[1]);
	MATFile *mfp = matOpen(obsFile.c_str(), "r");
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
		if (i%100 == 0) printf("Computing cooccurrences up to column %d...\n", i);
		for (j=i+1; j<c; j++) {
			Nuv[i][j] = 0;
			for (k=0; k<r; k++) Nuv[i][j] += (ip[i*r+k]&&ip[j*r+k]);
		}
	}

	string preFile(argv[2]); //file from buildMST.cpp
	mfp = matOpen(preFile.c_str(), "r");
	M = matGetVariable(mfp, "prev");

	int *prev;
	pr = mxGetPr(M);
	prev = (int *)mxMalloc(c * sizeof(*prev));
	for(i=0; i<c; i++) prev[i] = pr[i];

	//fabmap will have "run" trials and "tally" will keep track of image pairs
	int run = 100;
	int **tally;
	tally = new int*[r];
	for (i=0; i<r; i++) *(tally+i) = new int[r];
	for (i=0; i<r; i++) {
		for (j = i+1; j<r; j++) tally[i][j] = 0;
	}

	double *e, R, thresh = strtod(argv[3], NULL);
	R = r;

	e = new double[c]; //average of word presence
	for (i = 0; i < c; i++)		e[i] = N[i]/R;

	double detector[] = {0.88, 0.17, 0.12, 0.83};
	double *initialLocation;
	initialLocation = new double[c];
	for (i = 0; i < c; i++) initialLocation[i] = e[i];
	int offset = std::min<int>(r/6, 20)-1;

	printf("Loop closure trials to run: %d\n", run);
	for (int x = 0; x < run; x++) {
		printf("Trial %d: Computing marginal probabilities of the existence of words...\n", x+1);

		double **locations;
		locations = new double*[r];

		double *LP, num, tmp;
		int cur;
		LP = new double[r];
		for (i = 0; i < offset; i++) {
			locations[i] = new double[c];
			for (j = 0; j < c; j++) {
				tmp = initialLocation[j];
				cur = ip[j*r+i];
				num = detector[cur*2+1] * tmp;
				locations[i][j] = num/(num + detector[cur*2] * (1-tmp));
			}
		}

		int numLoop = 0, l, q, p, maxLoc, perc = 10, maxInd = -1, *loopArray;
		double bo, prob, aco, bco, beta, alpha, ratio, ratio0, ratio1, normal, maxProb;
		double *OL, *pdf;

		loopArray = new int[r-offset];

		printf("Analyzing images...");
		for (i = offset; i < r; i++) {
			if (i == 108) 
				printf("");
			int loopFlag = 0;

			if (100*i/R > perc) {
				printf("%2.f%% done...\n", 100*i/R);
				perc += 10;
			}

			cur = ip[i];
			const int s = i-offset-numLoop+1;
			OL = new double[s];
			maxProb = 0;
			for (l = 0; l < s; l++) {
				prob = detector[cur*2] * (1-locations[l][0]) + detector[cur*2+1] * locations[l][0];
				for (q = 1; q < c; q++) {
					p = prev[q];
					bo = (p > q) ? Nuv[q][p] : Nuv[p][q];
					/*if (p > q) {
						bo = Nuv[q][p];
					} else {
						bo = Nuv[p][q];
					}*/
					bo /= R;
					aco = (ip[p*r+i] == 0) ? (1-e[q]-e[p]+bo)/(1-e[p]) : (e[p]-bo)/e[p];
					bco = (ip[p*r+i] == 0) ? (e[q]-bo)/(1-e[p]) : bo/e[p];
					/*if (ip[p*r+i] == 0) {
						aco = (1-e[q]-e[p]+bo)/(1-e[p]);
						bco = (e[q]-bo)/(1-e[p]);
					} else {
						aco = (e[p]-bo)/e[p];
						bco = bo/e[p];
					}*/
					cur = ip[q*r+i];
					if (cur == 0) {
						//beta = e[q];
						tmp = aco;
						aco = bco;
						bco = tmp;
					} /*else {
						beta = 1-e[q];
					}*/
					beta = (cur == 0) ? e[q] : 1-e[q];
					alpha = 1-beta;

					ratio = alpha*aco/(beta*bco);
					ratio0 = detector[-2*cur+2]/detector[2*cur]*ratio;
					ratio1 = detector[-2*cur+3]/detector[2*cur+1]*ratio;

					prob *= (1-locations[l][q])/(1+ratio0) + locations[l][q]/(1+ratio1);
				}
				OL[l] = prob;
				if (prob > maxProb)	{
					maxProb = prob;
					maxInd = l;
				}
			}

			//sampling
			double M = std::max<int>(10, r/4);
			int *A = new int[M];
			fill_with_unique_rand(A, r, M);
			double total = 0, *samp;

			for (k = 0; k < M; k++) {
				int a = A[k];
				//avoid current image as well as image that likely forms a loop closure
				if (a != i && a != maxInd) {
					samp = new double[c];
					for (j = 0; j < c; j++) {
						tmp = initialLocation[j];
						cur = ip[j*r+a];
						num = detector[cur*2+1] * tmp;
						samp[j] = num/(num + detector[cur*2] * (1-tmp));
					}
					prob = detector[cur*2] * (1-samp[0]) + detector[cur*2+1] * samp[0];
					for (q = 1; q < c; q++) {
						p = prev[q];
						bo = (p > q) ? Nuv[q][p] : Nuv[p][q];
						bo /= R;
						aco = (ip[p*r+i] == 0) ? (1-e[q]-e[p]+bo)/(1-e[p]) : (e[p]-bo)/e[p];
						bco = (ip[p*r+i] == 0) ? (e[q]-bo)/(1-e[p]) : bo/e[p];
						cur = ip[q*r+i];
						if (cur == 0) {
							tmp = aco;
							aco = bco;
							bco = tmp;
						}
						beta = (cur == 0) ? e[q] : 1-e[q];
						alpha = 1-beta;

						ratio = alpha*aco/(beta*bco);
						ratio0 = detector[-2*cur+2]/detector[2*cur]*ratio;
						ratio1 = detector[-2*cur+3]/detector[2*cur+1]*ratio;

						prob *= (1-samp[q])/(1+ratio0) + samp[q]/(1+ratio1);
					}
					total += prob;
					delete[] samp;
				}
			}

			/*double avg = detector[cur*2] * (1-e[0]) + detector[cur*2+1] * e[0];
			for (q = 1; q < c; q++) {
				p = prev[q];
				bo = (p > q) ? Nuv[q][p] : Nuv[p][q];
				bo /= R;
				aco = (ip[p*r+i] == 0) ? (1-e[q]-e[p]+bo)/(1-e[p]) : (e[p]-bo)/e[p];
				bco = (ip[p*r+i] == 0) ? (e[q]-bo)/(1-e[p]) : bo/e[p];
				cur = ip[q*r+i];
				if (cur == 0) {
					tmp = aco;
					aco = bco;
					bco = tmp;
				}
				beta = (cur == 0) ? e[q] : 1-e[q];
				alpha = 1-beta;

				ratio = alpha*aco/(beta*bco);
				ratio0 = detector[-2*cur+2]/detector[2*cur]*ratio;
				ratio1 = detector[-2*cur+3]/detector[2*cur+1]*ratio;

				avg *= (1-e[q])/(1+ratio0) + e[q]/(1+ratio1);
			}*/

			for (l = 0, normal = 0; l < s; l++) normal += OL[l];
			normal = normal/s + 0.9 * total/M; //max(total/M, avg);

			pdf = new double[s];
			for (l = 0; l < s; l++) {
				pdf[l] = OL[l]/s/normal;
				if (l != s-1 && pdf[l] > thresh) { 
					loopFlag = 1;
					maxLoc = l;
				}
			}

			/*if (i == r-1) {
				for (int S = 0; S < M; S++) printf("%d ", A[S]);
				printf("\n");
				printf("total: %e\t avg: %e\t normal: %e\n", total, avg, normal);
			}*/
			if (loopFlag) {
				int count = 0;
				for (int x = 0; x < maxLoc - offset + numLoop; x++) {
					if (loopArray[x] == 1) count += 1;
				}
				printf("Loop closure detected with probability %f!\n", pdf[maxLoc]);
				printf("With observation %d (corresponds to image %d) and location %d (corresponds to image %d)\n", i, i+1, maxLoc, maxLoc+count+1);
				numLoop += 1;
				loopArray[i-offset] = 1;
				tally[maxLoc+count][i] += 1;
				//if (i == r-1)
				//	for (l = 0; l < s; l++) printf("%f\t", pdf[l]); //for(int S = 0; S < M; S++) printf("%d ", A[S]);
			} else {
				if (i == r-1) {
					for (l = 0; l < s; l++) printf("%f\t", pdf[l]);
					printf("\n");
				}
				loopArray[i-offset] = 0;
				locations[i-numLoop] = new double[c];
				for (j = 0; j < c; j++) {
					tmp = initialLocation[j];
					cur = ip[j*r+i];
					num = detector[cur*2+1] * tmp;
					locations[i-numLoop][j] = num/(num + detector[cur*2] * (1-tmp));
				}
			}

			delete[] pdf;
		}
		for (i=0; i<r-numLoop; i++) delete[] locations[i];
		delete[] locations;
	} 	
	
	priority_queue<Entry, vector<Entry>, CompareEntry> pq;

	for (i=0; i<r;i++) {
		for (j=i+1; j<r; j++) {
			if (tally[i][j] > 0) {
				Entry e1 = { j, i, tally[i][j] };
				pq.push(e1);
			}
		}
	}

	int pos = obsFile.find("obs");
	int posend = obsFile.find(".mat");
	string dataset = obsFile.substr(pos+3, posend-3);
	string result = "result";
	result.append(dataset);
	result += ".txt";

	FILE* file = NULL; 
	file = fopen(result.c_str(), "w");

	printf("Results written to %s\n", result.c_str());

	fprintf(file, "%d trials ran; probability threshold at %f\n", run, thresh);

	while(!pq.empty()) {
		Entry e1 = pq.top();
		pq.pop();
		fprintf(file, "%d\t-\t%d:\t%d\n", e1.obs+1, e1.loc+1, e1.count);
	}
	fclose(file);	

	return 0;

}

