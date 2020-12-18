#ifndef SETUP_H
#define SETUP_H
#include <vector>
#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/grid2d.h>

class Setup
{
public:
    Setup();
    Setup(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_);
    double calc_dt(double coef);
    std::vector<double> meshgridX();
    std::vector<double> meshgridY();
    std::vector<double> initial_condREN(std::vector<double> X, std::vector<double> Y);
    std::vector<double> initial_cond(std::vector<double> X, std::vector<double> Y); 
    std::vector<double> true_solution(std::vector<double> X, std::vector<double> Y, std::vector<double> VelX, std::vector<double> VelY, double dt);
    std::vector<double> true_solution2(std::vector<double> X, std::vector<double> Y, std::vector<double> VelX, std::vector<double> VelY, double dt);
    std::vector<double> initial_cond2(std::vector<double> X, std::vector<double> Y);
    std::vector<double> true_Laplacian(std::vector<double> X, std::vector<double> Y);

private:
    int N,M;
    double xmin,xmax,ymin,ymax;
};

#endif // SETUP_H
