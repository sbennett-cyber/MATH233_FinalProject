#ifndef NETWORKFINITEDIFFERENCE_H
#define NETWORKFINITEDIFFERENCE_H
#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/cf_2.h>
#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/sparsematrix_crs.h>
#include <iostream>
#include <vector>

class NetworkFiniteDifference
{
private:
    int NL;
    int NU;
    int Channels;
    double dx;
    double Start;
    double Mid;
    double End;
    int NNx;
    std::vector<double> L;
    std::vector<double> U;
    std::vector<double> x0;
    CF_2 *Initial;

public:
    void Setup(double dx_, int Channels_, double Start_, double Mid_, double End_);
    void GenerateMatrix(SparseMatrix_CRS &A, std::vector<double> &rhs, double r1, double r2, double r3);
    void InitialCondition(std::vector<double> &x);
    int return_NNx();
    void set_IC(CF_2 &IT);
};

#endif // NETWORKFINITEDIFFERENCE_H
