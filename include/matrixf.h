#include "in.h"
using namespace std;
#ifndef MATRIXF_H
#define MATRIXF_H

class Matrixf
{
public:
    void zMatrix(int nn, iMatrix& aa);
    void allMatrix(int nn,iMatrix& aa);
    void pMatrix(int nn, iMatrix& aa);
    void pMatrixD(int nn,dMatrix& aa);
    void CnMatrix(int nn, iMatrix& aa);
    void dD(int nn, iMatrix& aa);
    void vDeg(int nn, dRow &v, iMatrix &aa);
    double avrageVec(int nn, dRow &v);
    double varianceVec(int nn, dRow &v);
    void conMatrix(int nn,iMatrix &a,iMatrix &b);
    void DeleteMatElement(int nn,iMatrix &a);
    void mat_shrink_to_fit(int nn,iMatrix &a);
    void dD_eff(int nn, iMatrix& aa,int ans);
    void stepSumD(vector<double> &S, vector<double> &sum, int N, int stepL);

};

#endif // MATRIXF_H
