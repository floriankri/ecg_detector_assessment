//
//  minmaxfilter.cpp
//  Sortfilt
//
//  Created by Robert Weiss on 19/06/13.
//  Copyright (c) 2013 Robert Weiss. All rights reserved.
//

#include <algorithm>
#include <deque>
#include <float.h>
#include <iostream>
#include <vector>

#include "minmaxfilter.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) {
    double *inMatrix;
    int n, N;
    int rows, cols;
    bool max;

    inMatrix = mxGetPr(prhs[0]);
    n = mxGetScalar(prhs[1]);
    N = mxGetNumberOfElements(prhs[0]);

    max = mxGetScalar(prhs[2]) == 1;

    rows = (int) mxGetM(prhs[0]);
    cols = (int) mxGetN(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);

    minmaxfilter(mxGetPr(plhs[0]), inMatrix, n, N, max);
}

void minmaxfilter(double *result, double *x, int n, int N, bool max) {
    if (max)
        maxfilter(result, x, n, N);
    else
        minfilter(result, x, n, N);
}

void printDeque(std::deque<double> *D) {
    std::cout << "Deque: ";

    for (int i = 0; i < D->size(); i++) {
        std::cout << D->at(i) << " ";
    }

    std::cout << std::endl;
}

void minfilter(double *result, double* x, int n, int N) {
    int N1, N2, A, B, last_B, index;

    if (n % 2 == 0) {
		N1 = (n/2) - 1;
		N2 = (n/2);
	} else {
		N1 = (n-1) / 2;
		N2 = (n-1) / 2;
	}

    last_B = 0;
    index = -1;

    for (int i = 0; i < N; i++) {
        A = std::max(0, i-N1);
        B = std::min(N-1, i+N2);

        if (index < A) {
			double MIN = x[A];
			for (int j = A; j <= B; j++)
				if (x[j] <= MIN) {
					MIN = x[j];
					index = j;
				}
		}
		
		if (last_B != B)
			if (x[B] <= x[index])
				index = B;
		
		result[i] = x[index];
		
		last_B = B;
    }
}

void maxfilter(double *result, double* x, int n, int N) {
    int N1, N2, A, B, last_B, index;

    if (n % 2 == 0) {
		N1 = (n/2) - 1;
		N2 = (n/2);
	} else {
		N1 = (n-1) / 2;
		N2 = (n-1) / 2;
	}

    last_B = 0;
    index = -1;

    for (int i = 0; i < N; i++) {
        A = std::max(0, i-N1);
        B = std::min(N-1, i+N2);

        if (index < A) {
			double MAX = x[A];
			for (int j = A; j <= B; j++)
				if (x[j] >= MAX) {
					MAX = x[j];
					index = j;
				}
		}
		
		if (last_B != B)
			if (x[B] >= x[index])
				index = B;
		
		result[i] = x[index];
		
		last_B = B;
    }
}
