#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>


class Matrix
{
private:
public:
    Matrix() {};
    virtual void   mat_Vec_Product(const std::vector <double> &x, std::vector <double> &Ax) const = 0 ;
    virtual void   add_element(int i, int j, double v)        {} ;
    virtual double get_element(int i, int j)  const = 0 ;
    virtual ~Matrix() {};
};

#endif // MATRIX_H
