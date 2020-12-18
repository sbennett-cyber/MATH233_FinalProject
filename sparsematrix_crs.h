#ifndef SPARSEMATRIX_CRS_H
#define SPARSEMATRIX_CRS_H

#include <vector>
#include "matrix.h"

class SparseMatrix_CRS : public Matrix
{
private:
    std::vector <int> index;
    std::vector <double> values;
    std::vector <int> columns;
public:
    SparseMatrix_CRS();
    void mat_Vec_Product(const std::vector <double> &x, std::vector <double> &Ax) const;
    double get_element(int i, int j) const;
    void add_element(int i, int j, double v);
    int return_index(int i, int j) const;
    void print();
    ~SparseMatrix_CRS() {};
};

#endif // SPARSEMATRIX_CRS_H
