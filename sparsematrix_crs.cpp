#include "sparsematrix_crs.h"
#include <iostream>
#include <iomanip>


SparseMatrix_CRS::SparseMatrix_CRS()
{
    // We need this because no matter what the first val
    // of index should be one and if we are doiong
    // index.size - 1 the size can't be zero
    index.resize(1);
    index[0] = 0;
    values.resize(0);
    columns.resize(0);
}

void SparseMatrix_CRS::mat_Vec_Product(const std::vector <double> &x, std::vector <double> &Ax) const
{
    // get max column index
    int M = 0;
    for (int r = 0; r < columns.size(); ++r) {
        if (columns[r] > M)
        {
            M = columns[r];
        }
    }

    if(M >= x.size()){
        std::cout<<"Matrix and vector not same size"<<std::endl;
        return;
    }

    Ax.resize(index.size());
    for(int i = 0; i < index.size(); i++){
        int index_max = ( i == index.size()-1 ) ? columns.size()  : index[i+1]  ;
        Ax[i] = 0.;
        for(int j = index[i]; j < index_max; j++){
            Ax[i] = Ax[i] + values[j]*x[columns[j]];
//            std::cout<<values[j]<<std::endl;
//            std::cout<index[i]<<std::endl;
//            std::cout<<columns[j]<<std::endl;
        }
    }

}

int SparseMatrix_CRS::return_index(int i, int j) const
{
    // When i = 0, index size = 1
    if ( i >= index.size()  )
    {
       //Maxime @ Majerle: we don't need to print this. Also technically this is not an error, it just means that the corresponding coefficient is 0

//        std::cout << "ERROR: This is not in the matrix" << std::endl;
        return -1;
    }

    // If we are one from the max index the final value should be the column size.
    // We want to iterate through the columns of the last row which is
    // columns[index[i]] - columns[columns.size() - 1]

    //                  condition           ?  if true         :   if false
    int index_max = ( i == index.size()-1 ) ? columns.size()  : index[i+1]  ;

    // iteration through all columns in row i and do a column check
    for (int k = index[i]; k < index_max ; ++k) {
        if (columns[k] == j)
        {
            return k;
        }
    }
    return -1;
}

double SparseMatrix_CRS::get_element(int i, int j) const
{
    // get index
    int k = return_index(i, j);

    // we either return 0 or the value
    return ( k == -1) ? 0. : values[k] ;

}

void SparseMatrix_CRS::add_element(int i, int j, double v)
{
   // std::cout<<v<<std::endl;
    // get index
    int k = return_index(i, j);

    if (k != -1)
    {
        // if value already exists, add value to previous value
        values[k] += v;
    }
    else
    {
        if (i < index.size() - 1 )
        {
              // if you are not in the latest row, inform and do nothing
              std::cout << "ERROR: not the latest row" << std::endl;
        }
        if (i == index.size() - 1 )
        {
            // congrats, we're in the current row!
 //         std::cout << "current row " << std::endl;
            // push v to values and j to columns
            values.push_back(v);
            columns.push_back(j);
        }
        if (i > index.size() - 1 )
        {
            // we need to increase index!
            int current_size = index.size();
//            std::cout << "row +1 " << std::endl;

            // we iterate from current row to the row we want to be in!
            for (int l = current_size; l <= i; ++l)
            {
                // we push back columns.size to index
                // if we're adding a row with nothing we repeat same index.
                index.push_back( columns.size()  );
            }

            // when we get to correct row, push it push it real good.
            values.push_back(v);
            columns.push_back(j);



        }
    }

}


void SparseMatrix_CRS::print()
{
    // We don't need to +1 N, because 0 is the 1 element
    int N = index.size() ;

    int M = 0.;
    for (int r = 0; r < columns.size(); ++r) {
        if (columns[r] > M)
        {
            M = columns[r];
        }
    }
    // we need to +1 to M because we use <
    M += 1;

    std::cout << "N value " << N << std::endl;
    std::cout << "M value " << M << std::endl;

    for (int q = 0; q < N; ++q) {
        for (int p = 0; p < M; ++p) {
            std::cout << std::setprecision(3)<< get_element(q,p) << " ";
        }
         std::cout << std::endl;
    }

}
