
#include <iostream>
#include "networkfinitedifference.h"
#include <cmath>
#include <chrono>
#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/grid2d.h>
#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/matrix.h>
#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/sparsematrix_crs.h>
#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/conjugate_gradient.h>

using namespace std;

// Shayna Bennett
// Heat equation on a Y-Network

class velocity_X : public CF_2{
public:
    double operator()(double x, double y) const{
        return 2.*sin((3.1415*x)/2.)-sin(3.1415*x)+4*sin(2.*3.1415*x);
//        return 2*sin(x);
    }
};

NetworkFiniteDifference Net;


int main()
{
    double average=0.;
    double dx = 0.05;
    double dt = dx*dx;
    double n = 2.;
    double Start = 0;
    double Mid = 1;
    double End = 2;
//for (int Q = 0; Q<10; Q++)
//{

    vector<double> rhs;
    vector<double> sols;
    double D = 1./4.;
    velocity_X IC;
    SparseMatrix_CRS A;
    conjugate_gradient solver;


    Net.Setup(dx, 3, Start, Mid, End);
    Net.set_IC(IC);
    Net.InitialCondition(sols);
    int NNx = Net.return_NNx();

    rhs.resize(3*NNx);

    //Setup r=D/2h^2
    double r1 = (D*dt)/(dx*dx);
    double r2 = (D*dt)/(dx*dx);
    double r3 = (D*dt)/(dx*dx);

    //Create matrix of coefficents and initial guess
    Net.GenerateMatrix(A, rhs, r1, r2, r3);

    //Print initial condition
//    cout << "Solution at Initial Time" << endl;
//    for (int i = 0; i <3*NNx; i++)
//    {
//        cout<<sols[i]<<endl;
//    }


    //Solve in time
//    auto start = std::chrono::steady_clock::now();
    for (int i = 0; i <(n/dt); i++)
    {
        solver.solveBCG(A,sols, rhs);
    }
//    auto end = std::chrono::steady_clock::now();
//    std::chrono::duration<double> elapsed_seconds = end-start;
//    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";


   // A.print();
    cout << "Solution at Final Time" << endl;

    //Print Solution at Final Time
    for (int i = 0; i <3*NNx; i++)
    {
        cout<<sols[i]<<endl;
    }
//    average = average+elapsed_seconds.count();
//}
//    cout << "Average = "<< (average/10.)/(n/dt) << endl;
    cout << "Hello World!" << endl;
    return 0;
}
