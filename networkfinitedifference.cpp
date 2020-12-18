#include "networkfinitedifference.h"
#include <vector>
#include <iostream>

using namespace std;

void NetworkFiniteDifference::Setup(double dx_, int Channels_,double Start_, double Mid_, double End_)
{
 dx = dx_;
 Channels = Channels_;
 Start = Start_;
 Mid=Mid_;
 End= End_;
 NL=((Mid-Start)/dx)+1;
 NU=((End-Mid)/dx)+1;
 NNx=NL-2;

vector<double> dxLower;
vector<double> dxUpper;

 for (double i = Start; i < Mid; i+dx)
 {
     L.push_back(i);
     i=i+dx;
 }
 for (double i = Mid; i < End; i+dx)
 {
     U.push_back(i);
     i=i+dx;
 }

}

int NetworkFiniteDifference::return_NNx()
{
    return NNx;
}

void NetworkFiniteDifference::set_IC(CF_2 &IT){
    Initial = & IT;
}

void NetworkFiniteDifference::InitialCondition(vector<double> &x)
{
    for (double i = 1; i < NL-1; i++)
    {
        x.push_back((*Initial)(L[i],1));
    }
    for (double i = 1; i < NU-1; i++)
    {
        x.push_back((*Initial)(L[i],1));
    }
    for (double i = 1; i < NU-1; i++)
    {
        x.push_back((*Initial)(U[i],1));
    }

}

void NetworkFiniteDifference::GenerateMatrix(SparseMatrix_CRS &A, std::vector<double> &rhs, double r1, double r2, double r3)
{
    double cof = 2./9.;

    // Points to fix in junction are at NNx, NNx+1, and 2NNx+1

    // Iterate though all points on network
    for(int k = 0; k < 3*NNx; ++k)
    {
        if (k == 0 || k == NNx)                             // Boundary condition is zero, so remove left point
        {
            A.add_element(k,0,1.+2.*r1);
            A.add_element(k,1,-r1);
            rhs[k] = 1.;
        }
        if (k == 3*NNx-1 )      // Boundary condition is zero, so remove right point
        {
            A.add_element(k,k,1.+2.*r2);
            A.add_element(k,k-1,-r2);
            rhs[k] = 1.;
        }
        if (k == NNx-1)                         // Includes a point on the junction, so must be edited
        {
            A.add_element(k,NNx-2,-r1*(1.-0.5*cof));
            A.add_element(k,NNx-1,1.-r1*2.*(cof-1.));

            A.add_element(k,2*NNx-1,-r1*(2.*cof));
            A.add_element(k,2*NNx-2,r1*(0.5*cof));

            A.add_element(k,2*NNx,-r1*(2.*cof));
            A.add_element(k,2*NNx+1,r1*(0.5*cof));
            rhs[k] = 1.;
        }
        if (k == 2*NNx-1)                           // Includes a point on the junction, so must be edited
        {
            A.add_element(k,   NNx-2,   r2*((0.5*cof)));
            A.add_element(k,   NNx-1,   r2*((-2.*cof)));

            A.add_element(k, 2*NNx-1,1.-r2*2.*(cof-1.));
            A.add_element(k, 2*NNx-2,  -r2*(1.-0.5*cof));

            A.add_element(k, 2*NNx,     r2*(-2.*cof));
            A.add_element(k, 2*NNx+1,   r2*(0.5*cof));
            rhs[k] = 1.;
        }
        if (k == 2*NNx)                          // Includes a point on the junction, so must be edited
        {
            A.add_element(k,NNx-2,r3*((0.5*cof)));
            A.add_element(k,NNx-1,-r3*((2.*cof)));

            A.add_element(k,2*NNx,1.-r3*2.*(cof-1.));
            A.add_element(k,2*NNx+1,-r3*(1.-0.5*cof));

            A.add_element(k,2*NNx-1,r3*(-2*cof));
            A.add_element(k,2*NNx-2,r3*(0.5*cof));
            rhs[k] = 1.;
        }
        if (k > 0 && k<NNx-1)                    // Remaining points on Channel 1
        {
            A.add_element(k,k-1,-r1);
            A.add_element(k,k,1.+2.*r1);
            A.add_element(k,k+1,-r1);
            rhs[k] = 1.;
        }
        if (k > NNx && k<2*NNx-1)                // Remaining points on Channel 2
        {
            A.add_element(k,k-1,-r2);
            A.add_element(k,k,1.+2.*r2);
            A.add_element(k,k+1,-r2);
            rhs[k] = 1.;
        }
        if (k > 2*NNx && k<3*NNx-1)              // Remaining points on Channel 3
        {
            A.add_element(k,k-1,-r3);
            A.add_element(k,k,1.+2.*r3);
            A.add_element(k,k+1,-r3);
            rhs[k] = 1.;
        }
    }
}
