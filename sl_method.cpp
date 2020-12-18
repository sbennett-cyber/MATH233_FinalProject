#include "sl_method.h"
#include <iostream>
#include <cmath>
#include <omp.h>

SL_method::SL_method()
{

}

void SL_method::set_velX(CF_2 &velocity_x_t0){
    velocity_x = &velocity_x_t0;
}

void SL_method::set_velY(CF_2 &velocity_y_t0){
    velocity_y = &velocity_y_t0;
}

void SL_method::set_xd(std::vector<double> & xd_){
    xd=xd_;
}

void SL_method::set_yd(std::vector<double> & yd_){
    yd=yd_;
}

void SL_method::set_grid(Grid2D & grid_t0){
    grid = grid_t0;
}

void SL_method::set_solution(std::vector<double> & solution_t0_){
    solution = solution_t0_;
    solution_t0 = solution_t0_;

}

void SL_method::set_True(std::vector<double> & solution_true){
    True = solution_true;

}


void SL_method::Runge_Kutta2(double dt){
    int Size = grid.get_N()*grid.get_M();
    std::vector<double> Solution(Size);

#pragma omp parallel for
    for ( int n = 0; n<Size;n++ ){
        double x = grid.x_from_n(n);
        double y = grid.y_from_n(n);
        double x1 = x-(dt/2)*(*velocity_x)(x,y);
        double y1 = y-(dt/2)*(*velocity_y)(x,y);
        xd[n] = x - dt*(*velocity_x)(x1,y1);
        yd[n] = y - dt*(*velocity_y)(x1,y1);
    }
}

void SL_method::bilinear_interpolationENO(){

    std::vector<double> TempSol(grid.get_N()*grid.get_M());
    #pragma omp parallel for
    for ( int n = 0; n<grid.get_N()*grid.get_M();n++ ){

        //Get gridpoints to the right and left of x
        int i_min = (int) floor( (xd[n] - grid.get_xmin())/grid.get_dx() );
        i_min = std::max(0,i_min);
        int i_max = (int) ceil ((xd[n] - grid.get_xmin())/grid.get_dx());
        i_max = std::min(i_max, grid.get_N()-1);

        //Check if x is a gridpoint
        if (i_min == i_max){
            if (i_min==0){
                i_max = i_min+1;
            }
            else {
                i_min = i_max-1;
            }
        }

        //Get gridpoints above and below y
        int j_min = (int) floor( (yd[n] - grid.get_ymin())/grid.get_dy() );
        j_min = std::max(0,j_min);
        int j_max = (int) ceil ( (yd[n] - grid.get_ymin())/grid.get_dy() );
        j_max = std::min(j_max, grid.get_M()-1);

        //Check if y is a gridpoint
        if (j_min == j_max){
            if (j_min==0){
                j_max = j_min+1;
            }
            else {
                j_min = j_max-1;
            }
        }

        //Generate the location of the corners
        int Corner_00 = grid.n_from_ij(i_min, j_min);
        int Corner_01 = grid.n_from_ij(i_min, j_max);
        int Corner_10 = grid.n_from_ij(i_max, j_min);
        int Corner_11 = grid.n_from_ij(i_max, j_max);

        double x_min = grid.x_from_n(Corner_00);
        double y_min = grid.y_from_n(Corner_00);

        double x_max = grid.x_from_n(Corner_11);
        double y_max = grid.y_from_n(Corner_11);

        //Get dx and dy
        double dx = grid.get_dx();
        double dy = grid.get_dy();

        //Calculate the second derivatives at the corner points
        double DXX_00 = grid.dx_centered(solution, Corner_00);
        double DXX_01 = grid.dx_centered(solution, Corner_01);
        double DXX_10 = grid.dx_centered(solution, Corner_10);
        double DXX_11 = grid.dx_centered(solution, Corner_11);
        double DYY_00 = grid.dy_centered(solution, Corner_00);
        double DYY_01 = grid.dy_centered(solution, Corner_01);
        double DYY_10 = grid.dy_centered(solution, Corner_10);
        double DYY_11 = grid.dy_centered(solution, Corner_11);

        double MinX = grid.minmod(DXX_00,grid.minmod(DXX_01,grid.minmod(DXX_10,DXX_11)));
        double MinY = grid.minmod(DYY_00,grid.minmod(DYY_01,grid.minmod(DYY_10,DYY_11)));

        //Interpolate!
        TempSol[n] = (1./(dx*dy))*(solution[Corner_00]*(x_max-xd[n])*(y_max-yd[n]) +
                                   solution[Corner_01]*(x_max-xd[n])*(yd[n]-y_min) +
                                   solution[Corner_10]*(xd[n]-x_min)*(y_max-yd[n]) +
                                   solution[Corner_11]*(xd[n]-x_min)*(yd[n]-y_min))-
                                    ((xd[n]-x_min)*(x_max-xd[n])/2)*MinX-
                                    ((yd[n]-y_min)*(y_max-yd[n])/2)*MinY;

    }
    solution=TempSol;
    solution_t0 = TempSol;
}

void SL_method::reinitialize(double dt){
    // The method used here is from the paper "A second order accurate level set method on non-graded adaptive cartesian grids"
    std::vector<double> TempSol(grid.get_N()*grid.get_M());
    #pragma omp parallel for
    for ( int n = 0; n<grid.get_N()*grid.get_M();n++ )
    {

        //double S = (solution_t0[n]/(sqrt(solution_t0[n]*solution_t0[n] + grid.get_dx()*grid.get_dx())));
        double S = (solution_t0[n]/(sqrt(solution_t0[n]*solution_t0[n])));
        double DXB = grid.dx_backward(solution,n);
        double DXF = grid.dx_forward (solution,n);
        double DYB = grid.dy_backward(solution,n);
        double DYF = grid.dy_forward (solution,n);

        double H = Godunov_HT(solution_t0[n], DXF, DXB, DYF, DYB);

        TempSol[n] = solution[n] - dt*S*(H);

    }
    solution=TempSol;
}

double SL_method::Godunov_HT(double S, double dfx, double dbx, double dfy, double dby ){

    if (S<-0.0001)
    {
        double a = grid.max_double(dbx, 0.0);
        double b = grid.min_double(dfx, 0.0);
        double c = grid.max_double(dby, 0.0);
        double d = grid.min_double(dfy, 0.0);

        double H = sqrt(grid.max_double(a*a,b*b)+
                        grid.max_double(c*c,d*d))-1;
        return H;
    }
    else if (S>0.0001)
    {
        double a = grid.min_double(dbx, 0.0);
        double b = grid.max_double(dfx, 0.0);
        double c = grid.min_double(dby, 0.0);
        double d = grid.max_double(dfy, 0.0);

        double H = sqrt(grid.max_double(a*a,b*b)+
                        grid.max_double(c*c,d*d))-1;
        return H;
    }
    else
    {
        return 0.0;
    }
}

void SL_method::get_norms(){
    grid.two_Norm(True, solution);
    grid.Inf_Norm(True, solution);
    std::cout<<"Norms are calculated using ||True-Estimated||" <<std::endl;
}

void SL_method::save_vtk(std::string file_name){
    //    Grid2D grid(N,M, xmin, xmax, ymin, ymax);
    //save as .vtk file
    grid.initialize_VTK_file(file_name);
    grid.print_VTK_Format(True,"True",file_name);
    grid.print_VTK_Format(solution,"EstSol",file_name);

}
