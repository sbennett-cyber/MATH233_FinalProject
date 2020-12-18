#ifndef FVSOLVE_H
#define FVSOLVE_H

#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/grid2d.h>
#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/cf_2.h>
#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/sparsematrix_crs.h>
#include <vector>
#include <math.h>
class FVsolver
{
private:
    std::vector<double> solution;
    Grid2D grid;
    SparseMatrix_CRS matrix;
    std::vector<double> rhs;
    CF_2 *level_set;
    CF_2 *boundary_cond;
    CF_2 *forcing;
    double alpha;
    double mu;
    void build_linear_system();

public:
    FVsolver();
    inline void set_alpha(double alpha_in){alpha = alpha_in;}
    inline void set_mu(double mu_in){mu = mu_in;}
    inline void set_grid(Grid2D grid_in){grid = grid_in;}
    inline void set_level_set(CF_2 &ls_in){level_set = &ls_in;}
    inline void set_boundary_cond(CF_2 &bc_in){boundary_cond = &bc_in;}
    inline void set_forcing(CF_2 &force_in){forcing = &force_in;}
    std::vector<double> solve();
};

#endif // FVSOLVE_H
