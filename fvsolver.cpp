#include "fvsolver.h"

FVsolver::FVsolver()
{

}

void FVsolver::build_linear_system()
{
    rhs.resize(grid.get_N()*grid.get_M());

    for(int j = 0; j < grid.get_M(); j++){

        for(int i = 0; i < grid.get_N(); i++){
            double x_center = grid.x_from_n(grid.n_from_ij(i,j));
            double y_center = grid.y_from_n(grid.n_from_ij(i,j));

            // check the corners
            bool is_outside =
                    (*level_set)(x_center - 0.5*grid.get_dx(),y_center - 0.5*grid.get_dy()) > 0 &&
                    (*level_set)(x_center - 0.5*grid.get_dx(),y_center + 0.5*grid.get_dy()) > 0 &&
                    (*level_set)(x_center + 0.5*grid.get_dx(),y_center - 0.5*grid.get_dy()) > 0 &&
                    (*level_set)(x_center + 0.5*grid.get_dx(),y_center + 0.5*grid.get_dy()) > 0;

            if( is_outside ){
                matrix.add_element(grid.n_from_ij(i,j),grid.n_from_ij(i,j),1.);
                rhs[grid.n_from_ij(i,j)] = 0.;
            }
            else
            {
                // diagonal part of system first
                double size_x = (i != 0 && i != grid.get_N()-1) ? grid.get_dx() : 0.5*grid.get_dx();
                double size_y = (j != 0 && j != grid.get_M()-1) ? grid.get_dy() : 0.5*grid.get_dy();

                double volume_size = size_x * size_y;
                matrix.add_element(grid.n_from_ij(i,j),grid.n_from_ij(i,j), volume_size * alpha);
                rhs[grid.n_from_ij(i,j)] = volume_size * (*forcing)(x_center,y_center);

                // compute flux over the interface and add to rhs
                bool crossing =
                        (*level_set)(x_center - 0.5*grid.get_dx(),y_center - 0.5*grid.get_dy()) > 0 ||
                        (*level_set)(x_center - 0.5*grid.get_dx(),y_center + 0.5*grid.get_dy()) > 0 ||
                        (*level_set)(x_center + 0.5*grid.get_dx(),y_center - 0.5*grid.get_dy()) > 0 ||
                        (*level_set)(x_center + 0.5*grid.get_dx(),y_center + 0.5*grid.get_dy()) > 0;

                if (crossing){
                    rhs[grid.n_from_ij(i,j)] += sqrt(volume_size) * (*boundary_cond)(x_center,y_center);
                }

                double length_left = (
                            (*level_set)(x_center - 0.5*grid.get_dx(),y_center - 0.5*grid.get_dy()) > 0 &&
                            (*level_set)(x_center - 0.5*grid.get_dx(),y_center + 0.5*grid.get_dy()) > 0) ? 0 : grid.get_dy();

                double length_right = (
                            (*level_set)(x_center + 0.5*grid.get_dx(),y_center - 0.5*grid.get_dy()) > 0 &&
                            (*level_set)(x_center + 0.5*grid.get_dx(),y_center + 0.5*grid.get_dy()) > 0) ? 0 : grid.get_dy();

                double length_bottom = (
                            (*level_set)(x_center - 0.5*grid.get_dx(),y_center - 0.5*grid.get_dy()) > 0 &&
                            (*level_set)(x_center + 0.5*grid.get_dx(),y_center - 0.5*grid.get_dy()) > 0) ? 0 : grid.get_dx();

                double length_top = (
                            (*level_set)(x_center - 0.5*grid.get_dx(),y_center + 0.5*grid.get_dy()) > 0 &&
                            (*level_set)(x_center + 0.5*grid.get_dx(),y_center + 0.5*grid.get_dy()) > 0) ? 0 : grid.get_dx();

                if(i != 0){
                    matrix.add_element(grid.n_from_ij(i,j),grid.n_from_ij(i,j)  , mu * length_left / grid.get_dx());
                    matrix.add_element(grid.n_from_ij(i,j),grid.n_from_ij(i-1,j), -mu * length_left / grid.get_dx());
                }else{
                    rhs[grid.n_from_ij(i,j)] += length_left * (*boundary_cond)(x_center,y_center);
                }

                if(i != grid.get_N() - 1){
                    matrix.add_element(grid.n_from_ij(i,j),grid.n_from_ij(i,j)  , mu * length_right / grid.get_dx());
                    matrix.add_element(grid.n_from_ij(i,j),grid.n_from_ij(i+1,j), -mu * length_right / grid.get_dx());
                }else{
                    rhs[grid.n_from_ij(i,j)] += length_right * (*boundary_cond)(x_center,y_center);
                }

                if(j != 0){
                    matrix.add_element(grid.n_from_ij(i,j),grid.n_from_ij(i,j)  , mu * length_bottom / grid.get_dy());
                    matrix.add_element(grid.n_from_ij(i,j),grid.n_from_ij(i,j-1), -mu * length_bottom / grid.get_dy());
                }else{
                    rhs[grid.n_from_ij(i,j)] += length_bottom * (*boundary_cond)(x_center,y_center);
                }

                if(j != grid.get_M() - 1){
                    matrix.add_element(grid.n_from_ij(i,j),grid.n_from_ij(i,j)  , mu * length_top / grid.get_dy());
                    matrix.add_element(grid.n_from_ij(i,j),grid.n_from_ij(i,j+1), -mu * length_top / grid.get_dy());
                }else{
                    rhs[grid.n_from_ij(i,j)] += length_top * (*boundary_cond)(x_center,y_center);
                }

            }

        }
    }
}

std::vector<double> FVsolver::solve()
{
    build_linear_system();
    solution.resize(grid.get_N()*grid.get_M());

    matrix.print();

    return solution;
}
