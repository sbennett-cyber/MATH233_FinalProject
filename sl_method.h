#ifndef SL_METHOD_H
#define SL_METHOD_H
#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/grid2d.h>
#include <vector>
#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/cf_2.h>

class SL_method
{
private:
    Grid2D grid;
    std::vector<double> solution;
    std::vector<double> solution_t0;
    std::vector<double> xd;
    std::vector<double> yd;
    std::vector<double> True;
    CF_2 *velocity_x;
    CF_2 *velocity_y;


public:
    SL_method();

    double velocity_x_value(double x, double y);
    double velocity_y_value(double x, double y);
    double Godunov_HT(double S, double dfx, double dbx, double dfy, double dby );
    void save_vtk(std::string file_name);
    void bilinear_interpolationENO();
    void get_norms();
    void reinitialize(double dt);
    void set_xd(std::vector<double> & xd_);
    void set_yd(std::vector<double> & yd_);
    void Runge_Kutta2(double dt);
    void set_velX(CF_2 & velocity_x_t0);
    void set_velY(CF_2 & velocity_y_t0);
    void set_solution(std::vector<double> & solution_t0_);
    void set_True(std::vector<double> & solution_true);
    void set_grid(Grid2D & grid_t0);

};

#endif // SL_METHOD_H
