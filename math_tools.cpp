#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/math_tools.h>
#include <vector>
#include <cmath>

//math_tools::math_tools()
//{

//}

double math_tools::bilinear_interpolation(Grid2D & grid, const std::vector<double> & func, double x, double y){

    //Get gridpoints to the right and left of x
    int i_min = (int) floor( (x - grid.get_xmin())/grid.get_dx() );
        i_min = std::max(0,i_min);
    int i_max = (int) ceil ((x - grid.get_xmin())/grid.get_dx());
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
    int j_min = (int) floor( (y - grid.get_ymin())/grid.get_dy() );
        j_min = std::max(0,j_min);
    int j_max = (int) ceil ( (y - grid.get_ymin())/grid.get_dy() );
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



    //Interpolate!
    return (1./(dx*dy))*(func[Corner_00]*(x_max-x)*(y_max-y) +
                         func[Corner_01]*(x_max-x)*(y-y_min) +
                         func[Corner_10]*(x-x_min)*(y_max-y) +
                         func[Corner_11]*(x-x_min)*(y-y_min));


}


double math_tools::bilinear_interpolationENO(Grid2D & grid, const std::vector<double> & func, double x, double y){

    //Get gridpoints to the right and left of x
    int i_min = (int) floor( (x - grid.get_xmin())/grid.get_dx() );
        i_min = std::max(0,i_min);
    int i_max = (int) ceil ((x - grid.get_xmin())/grid.get_dx());
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
    int j_min = (int) floor( (y - grid.get_ymin())/grid.get_dy() );
        j_min = std::max(0,j_min);
    int j_max = (int) ceil ( (y - grid.get_ymin())/grid.get_dy() );
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



    //Interpolate!
    return (1./(dx*dy))*(func[Corner_00]*(x_max-x)*(y_max-y) +
                         func[Corner_01]*(x_max-x)*(y-y_min) +
                         func[Corner_10]*(x-x_min)*(y_max-y) +
                         func[Corner_11]*(x-x_min)*(y-y_min));


}

