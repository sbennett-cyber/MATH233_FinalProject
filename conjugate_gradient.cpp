#include "conjugate_gradient.h"
#include "matrix.h"
#include "sparsematrix_crs.h"
#include <vector>
#include <cmath>


conjugate_gradient::conjugate_gradient()
{

}

void conjugate_gradient::solveBCG(SparseMatrix_CRS & A, std::vector<double> & b, std::vector<double> & x)
{
    double TOL = 1E-10;
    double alpha = 1.;
    double omega = 1.;
    double resH = 1.;
    double res2 = 1.;
    double rhoOld = 1.;
    double rhoNew = 1.;
    double beta = 0.;

    std::vector<double> r;
    std::vector<double> p;
    std::vector<double> v;
    std::vector<double> h;
    std::vector<double> s;
    std::vector<double> t;
    std::vector<double> x0 = x;


    x.resize(b.size());
    p.resize(b.size());
    v.resize(b.size());
    h.resize(b.size());
    s.resize(b.size());
    t.resize(b.size());

    int size = x.size();

    // Calculate r=b-Ax
    A.mat_Vec_Product(x, r);
    for (int n=0; n < size; ++n)
    {
        r[n] = b[n] - r[n];
        p[n] = 0.;
        v[n] = 0.;
    }

    std::vector<double> rHat = r;

    for (int i = 0; i < size; i++)
    {
        rhoNew = scalarprod(rHat,r);
        beta = (rhoNew/rhoOld)*(alpha/omega);
        for (int n=0; n < size; ++n)
        {
            p[n] = r[n] + beta*(p[n]-omega*v[n]);
        }
        // Calculate v=Ap
        A.mat_Vec_Product(p, v);
        alpha = rhoNew/(scalarprod(rHat,v));
        for (int n=0; n < size; ++n)
        {
            h[n] = x[n] + alpha*p[n];
        }
        resH = scalarprod(r,r);
        if (std::sqrt(resH) < TOL)
        {
            x=h;
            break;
        }
        for (int n=0; n < size; ++n)
        {
            s[n] = r[n]- alpha*v[n];
        }
        A.mat_Vec_Product(s, t);
        omega = scalarprod(t,s)/scalarprod(t,t);
        for (int n=0; n < size; ++n)
        {
            x[n] = h[n] + omega*s[n];
        }
        res2 = scalarprod(s,s);
        if (std::sqrt(res2) < TOL)
        {
            break;
        }
        rhoOld=rhoNew;
        for (int n=0; n < size; ++n)
        {
            r[n] = s[n] - omega*t[n];
        }

    }
    b=x;
    x=x0;


}

void conjugate_gradient::solveCG(SparseMatrix_CRS & A, std::vector<double> & b, std::vector<double> & x)
{
    double TOL = 1E-10;
    double alpha = 0.;
    double res2 = 1.;
    std::vector<double> p;
    std::vector<double> Ap;
    double eps = 1E-20;

    x.resize(b.size());
    int size = x.size();

    std::vector<double> res ;
    std::vector<double> x0 = x;
    A.mat_Vec_Product(x, res);

    for (int n=0; n < size; ++n)
    {
        res[n] = b[n] - res[n];
    }

    p = res;

    while (std::sqrt(res2) > TOL)
    {
        A.mat_Vec_Product(p,Ap);
        res2 = scalarprod(res,res);
        alpha = res2/(std::max(scalarprod(Ap,p),eps));
        for (int n=0; n < size; ++n)
        {
            x[n] += alpha*p[n];
            res[n]  -= alpha*Ap[n];
        }
        double fac = scalarprod(res,res)/res2;
        for (int n=0; n < size; ++n)
        {
            p[n] = res[n] + fac*p[n];
        }
    }
    b=x;
    x0=x;
}

double conjugate_gradient::scalarprod(std::vector<double> &u, std::vector<double> &v)
{
    double product = 0;
    for (int k = 0; k < u.size(); ++k)
    {
        product = product+u[k]*v[k];
    }
    return product;
}
