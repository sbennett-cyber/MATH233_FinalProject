#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H
#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/grid2d.h>
#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/cf_2.h>
#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/fvsolver.h>


class conjugate_gradient
{
public:
    conjugate_gradient();
    void solveCG(SparseMatrix_CRS &A, std::vector<double> & b, std::vector<double> & x);
    void solveBCG(SparseMatrix_CRS &A, std::vector<double> & b, std::vector<double> & x);
    double scalarprod(std::vector<double> &u, std::vector<double> &v);
};

#endif // CONJUGATE_GRADIENT_H
