#include </Users/shaynabennett/MATH233_HOMEWORK/Important_Files/grid2d.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <omp.h>

using namespace std;

Grid2D::Grid2D()
{

}
Grid2D::Grid2D(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_)
{
    N=N_; //number of internal points
    M=M_; //number of internal points
    xmin=xmin_;
    xmax=xmax_;
    ymin=ymin_;
    ymax=ymax_;

    dx = (xmax - xmin) / (double) (N - 1);
    dy = (ymax - ymin) / (double) (M - 1);


}
double Grid2D::get_dx()
{
    return dx;
}

double Grid2D::get_dy()
{
    return dy;
}

double Grid2D::get_xmin()
{
    return xmin;
}

double Grid2D::get_ymin()
{
    return ymin;
}

int Grid2D::get_N()
{
    return N;
}

int Grid2D::get_M()
{
    return M;
}
int Grid2D::i_from_n(int n)
{
    //n = i + j * N
    return  n % M;
}
int Grid2D::j_from_n(int n)
{
    // Euclidean division, recall n = i + j * N
    return  n / M;
}
int Grid2D::n_from_ij(int i, int j)
{
    // Euclidean division, recall n = i + j * N
    int n = i+j*N;
    return  n;
}

double Grid2D::x_from_n(int n){
    return  xmin + dx*i_from_n(n);
}

double Grid2D::y_from_n(int n){
    return  ymin + dy*j_from_n(n);
}

double Grid2D::max_double(double a, double b){
    if (a>b){
        return a;
    }
    else{
        return b;
    }
}

double Grid2D::min_double(double a, double b){
    if (a<b){
        return a;
    }
    else{
        return b;
    }
}


double Grid2D::minmod(double D1, double D2){
    if ((D1>=0.0 && D2<=0.0) || (D1<=0.0 && D2>=0.0)){
        return 0.0;
    }
    else{
        if (abs(D1)<=abs(D2)){
            return D1;
        }
        else if (abs(D2)<=abs(D1)){
            return D2;
        }
        else{
            cout<<"Warning! Case not found in minmod. D1="<<D1<<" D2=" << D2 <<endl;
        }
    }
}

double Grid2D::minmod_Cell(double D1, double D2, double D3, double D4){
    if ((D1>=0.0 && D2<=0.0) || (D1<=0.0 && D2>=0.0) || (D3>=0.0 && D4<=0.0) || (D3<=0.0 && D4>=0.0)){
        return 0.0;
    }
    else{
        if (abs(D1)<=abs(D2)){
            return D1;
        }
        else if (abs(D2)<=abs(D1)){
            return D2;
        }
        else{
            cout<<"Warning! Case not found in minmod. D1="<<D1<<" D2=" << D2 <<" D3=" << D3<<" D4=" << D4<<endl;
        }
    }
}

//________________
double Grid2D::dx_forward(vector<double> &function, int n){
    try {
        // Check to see if gridpoint falls on the left boundary
        if (n % N <= 10e-7) {
            // First order discretization of the first derivative at the boundary and no flux BC, thus Qij = Qi-1j
            return 0.0;
        }
        else{

            double dxF = function.at(n)-function.at(n-1); //Solve Q_ij - Q_i-1j
            return dxF/dx;
        }
    }
    // Catch if value of n is outside the size of the vector function
    catch (const std::exception& e) {
        cout<<"Exception! n out of range."<< e.what() << '\n';
        return 0;
    }
}

\
double Grid2D::dx_backward(vector<double> &function, int n){
    try{
        //cout<<n-(N-1) << '\n';
        // Check to see if gridpoint falls on the right boundary
        if ((n+1) % N <= 10e-7){
            // First order discretization of the first derivative at the boundary and no flux BC, thus Qi+1j = Qij
            return 0.0;
        }
        else{
            double dxB = function.at(n+1)-function.at(n); //Solve Q_i+1j - Q_ij
            return dxB/dx;
        }
    }
    catch (const std::exception& e) {
        cout<<"Exception! n out of range."<< e.what() << '\n';
        return 0;
    }
}

double Grid2D::dy_forward(vector<double> &function, int n){
    try{
        if (n <= M-1){
            // First order discretization of the first derivative at the boundary and no flux BC, thus Qij = Qij-1
            return 0.0;
        }
        else{
            double dyF = function.at(n)-function.at(n-N); //Solve Q_ij - Q_ij-1
            return dyF/dy;
        }
    }
    catch (const std::exception& e) {
        cout<<"Exception! n out of range."<< e.what() << '\n';
        return 0;
    }
}

double Grid2D::dy_backward(vector<double> &function, int n){
    try{
        if (n <= M*N-1 && n>=(M*N)-(M)){
            // First order discretization of the first derivative at the boundary and no flux BC, thus Qij+1 = Qij
            return 0.0;
        }
        else{
            double dyB = function.at(n+N)-function.at(n); //Solve Q_ij+1 - Q_ij
            return dyB/dy;
        }


    }  catch (const std::exception& e) {
        cout<<"Exception! n out of range."<< e.what() << '\n';
        return 0;

    }
}

double Grid2D::dy_centered(vector<double> &function, int n){
    try{
        if (n <= M-1 && n>=0){
            // First order discretization of the first derivative at the boundary and no flux BC, thus Qij+1 = Qij
            double dyB = -function.at(n)+function.at(n+M); //Solve Q_ij+1 - Q_ij
            return dyB/(dy*dy);
        }
        else if (n <= M*N-1 && n>=(M*N)-(M)){
            // First order discretization of the first derivative at the boundary and no flux BC, thus Qij+1 = Qij
            double dyB = -function.at(n)+function.at(n-M); //Solve Q_ij+1 - Q_ij
            return dyB/(dy*dy);
        }
        else if (n > M*N-1 || n<0){
            // First order discretization of the first derivative at the boundary and no flux BC, thus Qij+1 = Qij
            return 0.0;
        }
        else{
            double dyB = function.at(n-M)-2*function.at(n)+function.at(n+M); //Solve Q_ij+1 - Q_ij
            return dyB/(dy*dy);
        }


    }  catch (const std::exception& e) {
        cout<<"Exception! n out of range in Y."<< e.what() << '\n';
        return 0;

    }
}

double Grid2D::dx_centered(vector<double> &function, int n){
    try{
        if (n % N <= 10e-7){
            // First order discretization of the first derivative at the boundary and no flux BC, thus Qij+1 = Qij
            double dyB = -function.at(n)+function.at(n+1); //Solve Q_ij+1 - Q_ij
            return dyB/(dy*dy);
        }
        else if ((n+1) % N <= 10e-7){
            // First order discretization of the first derivative at the boundary and no flux BC, thus Qij+1 = Qij
            double dyB = -function.at(n)+function.at(n-1); //Solve Q_ij+1 - Q_ij
            return dyB/(dy*dy);
        }
        else{
            double dyB = function.at(n-1)-2*function.at(n)+function.at(n+1); //Solve Q_ij+1 - Q_ij
            return dyB/(dy*dy);
        }


    }  catch (const std::exception& e) {
        cout<<"Exception! n out of range in X."<< e.what() << '\n';
        return 0;

    }
}


// initialize the .vtk file at the specified address with all the grid information
void Grid2D::initialize_VTK_file(std::string file_name){
    int node_of_cell[4];
    FILE *outFile = fopen(file_name.c_str(),"w");
    fprintf(outFile,"# vtk DataFile Version 2.0 \n");
    fprintf(outFile,"Quadtree Mesh \n");
    fprintf(outFile,"ASCII \n");
    fprintf(outFile,"DATASET UNSTRUCTURED_GRID \n");

    //% first output the list of nodes
    fprintf(outFile,"POINTS %d double \n",N*M);
    for (int n=0; n<N*M; n++){
        fprintf(outFile,"%e %e %e\n",x_from_n(n), y_from_n(n), 0.0);
    }
     // then output the list of cells. each cell is composed of four nodes
     fprintf(outFile,"CELLS %d %d \n",(N-1)*(M-1),5*(N-1)*(M-1));
     for (int i=0; i<N-1; i++){
         for (int j=0; j<M-1; j++) {
             node_of_cell[0] = n_from_ij(i  ,j  );
             node_of_cell[1] = n_from_ij(i+1,j  );
             node_of_cell[2] = n_from_ij(i+1,j+1);
             node_of_cell[3] = n_from_ij(i  ,j+1);
             fprintf(outFile,"%d %d %d %d %d\n",4,node_of_cell[0], node_of_cell[1], node_of_cell[2], node_of_cell[3]);
         }
     }
     fprintf(outFile,"CELL_TYPES %d \n",(N-1)*(M-1));

    for (int n=0; n<(N-1)*(M-1); n++)
    {
        fprintf(outFile,"%d \n",9);
    }
    fprintf(outFile,"POINT_DATA %d \n",N*M);
    fclose (outFile);

}

// this function write the values of the vector F into the vtk file. Before using it, the .vtk file must have been initialized with all the grid infos
void Grid2D::print_VTK_Format( std::vector<double> &F, std::string data_name, std::string file_name ){
    FILE *outFile;
    outFile = fopen(file_name.c_str(),"a");
    fprintf(outFile,"SCALARS %s double 1 \n",data_name.c_str());
    fprintf(outFile,"LOOKUP_TABLE default \n");
    for (int n=0; n<N*M; n++) {
        fprintf(outFile,"%e \n",F[n]);
    }
    fclose (outFile);
}

// Infinity Norm
void Grid2D::Inf_Norm(vector<double> True, vector<double> Est){
    vector<double> Diff;
    for (int i = 0; i<=True.size(); i++){
        Diff.push_back( True[i]-Est[i] );
    }
    cout << "The inf norm is " << *std::max_element(Diff.begin(),Diff.end()) << endl;
}

// 2 Norm
void Grid2D::two_Norm(vector<double> True, vector<double> Est){
    double Diff;
    for (int i = 0; i<=True.size(); i++){
        Diff=Diff+abs(True[i]-Est[i])*abs(True[i]-Est[i]);
    }
    cout << "The two norm is " << dx*dy*sqrt(Diff) << endl;
}

void Grid2D::print(){
    cout<<"N"<<N<< '\n';
    cout<<"M"<<M << '\n';
    cout<<"xmin"<<xmin << '\n';
    cout<<"xmax"<<xmax<< '\n';
    cout<<"ymin"<<ymin << '\n';
    cout<<"ymax"<<ymax<< '\n';

}

