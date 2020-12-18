#ifndef MATH_TOOLS_H
#define MATH_TOOLS_H
#include <vector>
#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/grid2d.h>

class math_tools
{
public:
   // math_tools();
    double bilinear_interpolation(Grid2D & grid, const std::vector<double> & func, double x, double y);
    double bilinear_interpolationENO(Grid2D & grid, const std::vector<double> & func, double x, double y);
};

#endif // MATH_TOOLS_H
