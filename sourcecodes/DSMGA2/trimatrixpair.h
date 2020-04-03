#ifndef _TRI_MATRIX_PAIR_
#define _TRI_MATRIX_PAIR_

#include <cassert>
#include <cstdio>
#include <utility>

using namespace std;

class TriMatrixPair {

public:
	TriMatrixPair() {
        matrix = NULL;
        size = 0;
    }
	TriMatrixPair(int n) {
        matrix = NULL;
        size = 0;
        init(n);
    }
    void init(int n) {
        if (n == 0)
            return;

        if (n == size) {
            for (int i = 0; i < size - 1; i++)
                for (int j = 0; j <= i; j++)
                    matrix[i][j] = pair<double, double>(0.0, 0.0);
                    //matrix[i][j] = T(0.0);
        } else {
            release();
            size = n;
            matrix = new pair<double, double>*[size - 1];
            for (int i = 0; i < size - 1; i++)
                matrix[i] = new pair<double, double>[i + 1];
            for (int i = 0; i < size - 1; i++)
                for (int j = 0; j <= i; j++)
                    matrix[i][j] = pair<double, double>(0.0, 0.0);
                    //matrix[i][j] = T(0.0);
        }
    }
    void release() {
        if (matrix == NULL)
            return;
        for (int i = 0; i < size - 1; i++)
            delete[] (matrix[i]);
        delete[] matrix;
        size = 0;
    }

    ~TriMatrixPair() {
        release();
    }
/*
    void dec(int i, int j, int step=1) {
        if (i==j) return;
        if (i < j)
            matrix[j-1][i]-=step;
        else
            matrix[i-1][j]-=step;
    }

    void inc(int i, int j, int step=1) {
        if (i==j) return;
        dec(i, j, -step);
    }
*/

    //void write(int i, int j, T val) {
    void write(int i, int j, pair<double,double> val) {
        assert(i < size && j < size);
        if (i == j) return;
        if (i < j) {
            int temp = i;
            i = j;
            j = temp;
        }
        matrix[i - 1][j] = val;
    }

    //T operator()(int i, int j) const {
    pair<double, double> operator()(int i, int j) const {
        if (i == j)
            return pair<double, double>(1.0, 1.0);
            //return T(1.0);
        assert(i < size && j < size);
        if (i < j) {
            int temp = i;
            i = j;
            j = temp;
        }
        return matrix[i - 1][j];
    }

private:
    pair<double, double>** matrix;
    //T** matrix;
	int size;

};
#endif
