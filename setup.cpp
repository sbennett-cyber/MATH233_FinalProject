#include "setup.h"
#include <iostream>
#include <cmath>
#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/grid2d.h>
#include <vector>

using namespace std;
Setup::Setup()
{

}

Setup::Setup(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_)
{

    N=N_; //number of internal points
    M=M_; //number of internal points
    xmin=xmin_;
    xmax=xmax_;
    ymin=ymin_;
    ymax=ymax_;
    //N_in=N-2;
    //M_in=M-2;

}

// Meshgrid in X //
vector<double> Setup::meshgridX(){
    Grid2D grid(N,M, xmin, xmax, ymin, ymax);
    double dx = grid.get_dx();
    vector<double> X;
    double valX = xmin+dx;
    //cout<<valX<<endl;

    for (int j=0;j<M;j++){
        for (int i=0;i<N;i++){
            X.push_back(valX);
            //cout<<valX<<endl;
            valX=valX+dx;
        }
        valX = xmin+dx;
    }
    return X;
}
//______________//
// Meshgrid in Y
vector<double> Setup::meshgridY(){
    Grid2D grid(N,M, xmin, xmax, ymin, ymax);
    double dx = grid.get_dx();
    vector<double> Y;
    double valY = ymin+dx;

    for (int j=0;j<M;j++){
        for (int i=0;i<N;i++){
            Y.push_back(valY);
            //cout<<valY<<endl;
        }
    valY=valY+dx;
    }
    return Y;
}
//__________________//
// Initial Condition
vector<double> Setup::initial_cond(vector<double> X, vector<double> Y){
    vector<double> t0;
    for (int j=0;j<N*M;j++){
        if ((sqrt((X[j]-0.5)*(X[j]-0.5)+Y[j]*Y[j]))-0.2<=0){
            t0.push_back(1);
        }
        else{
            t0.push_back(0);
        }

    }
     return t0;
}

//__________________//
// Initial Condition 2
vector<double> Setup::initial_cond2(vector<double> X, vector<double> Y){
    vector<double> t0;
    for (int j=0;j<N*M;j++){
        t0.push_back(((sqrt((X[j]-0.25)*(X[j]-0.25)+Y[j]*Y[j]))-0.2));
    }
     return t0;
}

// Initial Condition 2
vector<double> Setup::initial_condREN(vector<double> X, vector<double> Y){
    vector<double> t0;
    for (int j=0;j<N*M;j++){
        //t0.push_back(49*((sqrt((X[j]-0.25)*(X[j]-0.25)+Y[j]*Y[j]))-0.2));
        t0.push_back((X[j]*X[j]+Y[j]*Y[j])*(X[j]*X[j]+Y[j]*Y[j])-X[j]*X[j]+Y[j]*Y[j]);
    }
     return t0;
}

//______________//
// True Solution
vector<double> Setup::true_solution(vector<double> X, vector<double> Y, vector<double> VelX, vector<double> VelY, double dt){
    vector<double> t0;
    for (int j=0;j<N*M;j++){
        if ((sqrt((X[j]-0.5)*(X[j]-0.5)+Y[j]*Y[j]))-0.2<=0){
            t0.push_back(1);
        }
        else{
            t0.push_back(0);
        }

    }
     return t0;
}

// True Solution
vector<double> Setup::true_solution2(vector<double> X, vector<double> Y, vector<double> VelX, vector<double> VelY, double dt){
    vector<double> t0;
    for (int j=0;j<N*M;j++){
        t0.push_back(((sqrt((X[j]-0.25)*(X[j]-0.25)+Y[j]*Y[j]))-0.2));
    }
     return t0;
}

// True Solution
vector<double> Setup::true_Laplacian(vector<double> X, vector<double> Y){
    vector<double> t0;
    for (int j=0;j<N*M;j++){
        double psi_X = (X[j]-0.25) / sqrt((X[j]-0.25)*(X[j]-0.25)+Y[j]*Y[j]);
        double psi_Y = Y[j]/ sqrt((X[j]-0.25)*(X[j]-0.25)+Y[j]*Y[j]);
        t0.push_back(-Y[j]*psi_X+X[j]*psi_Y);
    }
     return t0;
}

//________________________//
// Calculate and return dt
double Setup::calc_dt(double coef){
    Grid2D grid(N,M, xmin, xmax, ymin, ymax);
    double dx = grid.get_dx();
    cout << "The value for dx = " << dx << " and dt =  " << coef*dx <<endl;
    return coef*dx;
}

