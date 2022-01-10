#include <iostream>
#include "CSRMat.h"
#include <memory>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <iterator>
#include <cmath>

using namespace  std;


template <class T>
CSRMat<T>::CSRMat(int rows, int cols, int nnzs, bool usr_provided_input_output)
{
    // set the values of this instance of the class as the values defined in the header
    this->rows = rows;
    this->cols = cols;
    this->no_values = rows*cols;
    this->nnzs = nnzs;
    this->usr_provided_input_output = usr_provided_input_output;

    // if (!this->usr_provided_input_output)
    // {
    //     // This is the condition to instantiate vectors if the user does not provide a 
    //     // input to the class and wants an automatically generated random CSR matrix 
    //     // Works in combination with the "Populate" function, but not yet completed
    // }
}

template <class T>
CSRMat<T>::CSRMat(int rows, int cols, int nnzs, std::vector<T>* values_ptr, std::vector<int>* row_position_ptr, std::vector<int>* col_index_ptr, std::vector<T>* b_ptr, std::vector<T>* output_sparse_ptr)
{
    // This is the overloaded constructor for this class and allows the user to directly instantiate
    // the class with their own vectors using pointers
    // set the values of this instance of the class as the values defined in the header
    this->rows = rows;
    this->cols = cols;
    this->no_values = rows*cols;
    this->nnzs = nnzs;
    this->usr_provided_input_output = usr_provided_input_output;

    // set the CSR pointers instantiated in the class as the ones given by the user
    this->values = values_ptr;
    this->row_position = row_position_ptr;
    this->col_index = col_index_ptr;
    this->b = b_ptr;
    this->output_sparse = output_sparse_ptr;
}

template <class T>
CSRMat<T>::~CSRMat()
{
    // as we are not creating new memory there is no need for special delete dunctions
    // in the destructor of the class
}

template <class T>
void CSRMat<T>::PrintMatrix_A()
{
    // Printing in CSR format which outputs the values, row_position and col_index vectors to the 
    // standard error channel on the console
    std::cerr << "Printing matrix A" << std::endl;
    std::cerr << "Values: ";
    for (int j=0; j<this->nnzs; j++)
    {
        std::cerr << this->values->at(j) << ", ";
    }
    std::cerr << std::endl;
    std::cerr << "row position: ";
    for (int j=0; j<this->rows + 1; j++)
    {
        std::cerr << this->row_position->at(j) << ", ";
    }
    std::cerr << std::endl;
    std::cerr << "column index: ";
    for (int j=0; j<this->nnzs; j++)
    {
        std::cerr << this->col_index->at(j) << ", ";
    }
    std::cerr << std::endl;
}

template <class T>
void CSRMat<T>::PrintVector_b()
{
    std::cerr << "Printing Vector b: " << std::endl;
    for (int i=0; i < this->rows; i++)
    {
        std::cerr << this->b->at(i) << ", ";
    }
    std::cerr << std::endl << std::endl;
}

template <class T>
void CSRMat<T>::PrintVector_x()
{
    std::cerr << "Printing Vector x: " << std::endl;
    for (int i=0; i < this->rows; i++)
    {
        std::cerr << this->output_sparse->at(i) << ", ";
    }
    std::cerr << std::endl << std::endl;
}

template <class T>
void CSRMat<T>::sparse2dense(T *dense_mat_ptr)
{
    // iterate over the row position
    for (int r=0; r<this->rows ; r++)
    {
        // iterate over the values whose index's are between the current pairing in the row_position
        for (int i=this->row_position->at(r) ; i<this->row_position->at(r+1); i++)
        {
            dense_mat_ptr[(r) * this->cols + this->col_index->at(i)] = this->values->at(i);
        }
    }
}

template <class T>
void CSRMat<T>::jacobi(double tolerance)
{   
    // fill the output vector with 0s to initialise it (wipe any previous values)
    std::fill(this->output_sparse->begin(), this->output_sparse->end(), 0);
    // check matrix A is square
    if (this->rows != this->cols)
    {
        std::cerr << "Matrix A is not square, pass in a square matrix" << std::endl;
        return;
    }

    // note the elements in a vector constructed this way are stored on the heap and the
    // header info to access it is on the stack
    // This will automatically be deleted once it goes out of scope at the end of the method
    std::vector<T> matrix_x1;
    matrix_x1.reserve(rows);

    int iterations = 0;
    int flag = 0;
    while (iterations < this->rows)
    {
        // iterate over the number of rows
        for (int r=0; r< this->rows; r++)
        {
            T sum_holder = 0;
            // for each row calculate what location the the values for this row are stored
            // I.e. the beginning and ending index of the row in the col_idx & value vector
            int pair_start = this->row_position->at(r);
            int pair_end = this->row_position->at(r+1);
            // find the value of the diagonal element in this row
            auto diagonal_location = std::find(this->col_index->begin() + pair_start, this->col_index->begin() + pair_end, r);
            int diagonal_idx = std::distance(this->col_index->begin(), diagonal_location);
            
            // this is the same as the normal jacobi algorithm but instead of iterating
            // over all values in the current row we just use the non zero values in the current row
            // which are located between (pair_start & pair_end) in the col_idx and values vector
            for  (int j = pair_start; j<pair_end; j++)
            {
                // skip the sum_holder if on the diagonal
                if (j == diagonal_idx) continue;
                {
                    // perform the summation part of Jacobi's method using the old guess (matrix_x1)
                    // element a_{i,j} * x_old{j}
                    sum_holder += this->values->at(j) * matrix_x1[this->col_index->at(j)];
                }

            }
            // update answer with (1 / diagonal element) * (b - sum_holder)
            this->output_sparse->at(r) = (1 / this->values->at(diagonal_idx)) *  (this->b->at(r) - sum_holder);
            // compare with previous value and check if below the tolerance
            if (abs(matrix_x1[r] - this->output_sparse->at(r)) < tolerance)
                flag++;
            matrix_x1[r] = this->output_sparse->at(r);
        }
        iterations ++;
    }
    std::cerr << "Iterations: " << iterations;
    
}

template <class T>
void CSRMat<T>::successive_over_relaxation(double tolerance, double r_factor)
{
    // fill the output vector with 0s to initialise it (wipe any previous values)
    std::fill(this->output_sparse->begin(), this->output_sparse->end(), 0);
    // this->output_sparse.reset(new double[this->rows] {0});
    // check matrix A is square
    if (this->rows != this->cols)
    {
        std::cerr << "Matrix A is not square, pass in a square matrix" << std::endl;
        return;
    }
    // will be deleted automatically once out of scope
    std::vector<T> matrix_x1;
    matrix_x1.reserve(rows);

    int iterations = 0;
    int flag = 0;
    while (iterations < this->rows)
    {
        // iterate over the number of rows
        for (int r=0; r< this->rows; r++)
        {
            T sum_holder = 0;
            // for each row calculate what location the the values for this row are stored
            // I.e. the beginning and ending index of the row in the col_idx & value vector
            int pair_start = this->row_position->at(r);
            int pair_end = this->row_position->at(r+1);
            // find the value of the diagonal element in this row
            auto diagonal_location = std::find(this->col_index->begin() + pair_start, this->col_index->end() + pair_end, r);
            int diagonal_idx = std::distance(this->col_index->begin(), diagonal_location);
            
            // this is the same as the normal jacobi algorithm but instead of iterating
            // over all values in the current row we just use the non zero values in the current row
            // which are located between (pair_start & pair_end) in the col_idx and values vector
            for  (int j = pair_start; j<pair_end; j++)
            {
                // skip the sum_holder operation if on the diagonal element
                if (j == diagonal_idx) continue;
                {
                    sum_holder += this->values->at(j) * matrix_x1[col_index->at(j)];
                }
            }
            // update answer with SOR formula
            this->output_sparse->at(r) += r_factor *((this->b->at(r) - sum_holder) / this->values->at(diagonal_idx) - this->output_sparse->at(r));
            // compare with previous value and check if below the tolerance
            if (abs(matrix_x1[r] - this->output_sparse->at(r)) < tolerance)
                flag++;
            matrix_x1[r] = this->output_sparse->at(r);
        }
        iterations ++;
    }
    std::cerr << "Iterations: " << iterations;
}

template <class T>
void CSRMat<T>::gauss_seidel(double tolerance)
{
    std::fill(this->output_sparse->begin(), this->output_sparse->end(), 0);
    // this->output_sparse.reset(new double[rows] {0});

    std::vector<T> matrix_x1;
    matrix_x1.reserve(rows);

    int iterations = 0;
    int flag = 0;
    while (flag<this->rows)
    {
        // iterate over the number of rows
        for (int r=0; r<this->rows; r++)
        {
            // for each row calculate what location the the values for this row are stored
            // I.e. the beginning and ending index of the row in the col_idx & value vector
            int pair_start = this->row_position->at(r);
            int pair_end = this->row_position->at(r+1);

            // find the value of the diagonal element in this row
            auto diagonal_location = std::find(this->col_index->begin() + pair_start, this->col_index->end() + pair_end, r);
            int diagonal_idx = std::distance(this->col_index->begin(), diagonal_location);

            T copy_current_guess = this->output_sparse->at(r);
            // rth element in b divided by the diagonal in A
            this->output_sparse->at(r) = this->b->at(r) / this->values->at(diagonal_idx);

            // this is the same as the normal G-S algorithm but instead of iterating
            // over all values in the current row we just use the non zero values in the current row
            // which are located between (pair_start & pair_end) in the col_idx and values vector
            for  (int j = pair_start; j<pair_end; j++)
            {
                if (j == diagonal_idx)
                    continue;

                this->output_sparse->at(r) = this->output_sparse->at(r) - ((this->values->at(j) / this->values->at(diagonal_idx)) * matrix_x1[this->col_index->at(j)]);
                matrix_x1[r] = this->output_sparse->at(r);
               
            }
            if (abs(copy_current_guess - this->output_sparse->at(r)) < tolerance)
                flag +=1;
        }
        iterations ++;
    }
    std::cerr << "Iterations: " << iterations;

}

template <class T>
void CSRMat<T>::Conjugate_Grad(double tolerance)
{
    // makse sure the output is set to 0
    std::fill(this->output_sparse->begin(), this->output_sparse->end(), 0);

    // Store them as pointers so that they can be passed easily to the matrix vector multiplier function
    // This actually stores both the header and contents of the vector on the heap but both will
    // be automatically destroyed when the vector goes out of scope delete the pointer to it at the end!
    std::vector<T>* r = new std::vector<T>(rows, 0); // residual vector
    std::vector<T>* d = new std::vector<T>(rows, 0); // the update to our guess
    std::vector<T>* q = new std::vector<T>(rows, 0); // this will be the product of A . d
    std::vector<T>* temp2 = new std::vector<T>(rows, 0); // this will be product of A . x

    // auto *r = new double[this->rows] {0}; // residual vector
    // auto *d = new double[this->rows] {0}; // the update to our guess
    // auto *q = new double[this->rows] {0}; // this will be the product of A . d
    // auto *temp2 = new double[this->rows] {0}; // this will be product of A . x
    T sigma_new = 0; // this is the value of the residual
    T sigma_0 = 0;
    T alpha = 0; // this is the length of the step we will take
    T sigma_old = 0; // this is to store the original sigma at each loop
    T beta = 0; // this is the ratio of the new residual : old residual

    // initialise "r" following r = b - Ax
    matA_vec_mult(this->output_sparse, r);; // r = Ax
    for (int i=0; i<this->rows; i++)
    {
        r->at(i) = this->b->at(i) - r->at(i);
    }
    //initialise d by setting it to r
    for (int i=0; i<this->rows; i++)
    {
        d->at(i) = r->at(i);
    }
    // initialise sigma_new by setting it equal to r^T . r
    sigma_new = vec_Trans_vec_mult(r, r);

    // initialise sigma_0 with the value of sigma_new
    sigma_0 = sigma_new;

    int iterations = 0;
    while ((iterations <  5000000) && (sigma_new > pow(tolerance, 2) * sigma_0))
    {
        // q = A . d
        matA_vec_mult(d, q);

        // calculate d^T . q and assign it to temp
        T temp = 0;
        temp = vec_Trans_vec_mult(d, q);
        // initialise alpha: alpha = sigma_new /  temp
        alpha = sigma_new / temp;

        // update x with our new approximation
        // x = x + alpha * d
        for (int i=0; i<this->rows; i++)
        {
            this->output_sparse->at(i) += alpha * d->at(i);
        }

        // update the residual
        // r = b - Ax
        // calculate Ax and assign it to temp2
        matA_vec_mult(this->output_sparse, temp2);
        for (int i=0; i<this->rows; i++)
        {
            r->at(i) = this->b->at(i) - temp2->at(i);
        }

        // update sigma_old
        sigma_old = sigma_new;

        // update sigma new = r^T . r
        sigma_new = vec_Trans_vec_mult(r, r);

        // calculate beta = sigma_new / sigma_old
        beta = sigma_new / sigma_old;

        // update d = r + beta * d
        for (int i=0; i< this->rows; i++)
        {
            d->at(i) = r->at(i) + beta * d->at(i);
        }
        iterations ++;
    }
    std::cerr << "Iterations: " << iterations;
    delete r;
    delete d;
    delete q;
    delete temp2;
    
}

template <class T>
void CSRMat<T>::Cholesky()
{
    // the structure of the cholesky matrix (L) is not always the exact same as the the original 
    // Matrix lower elements, Since we are stoting the matrix sparsely this means we first need to
    // determine where there will be elements in the cholesky matrix
    // This method can only solve Symmetric Matrices!

    // Reference: Solution of sparse positive definite systems on a hypercube
    // https://core.ac.uk/download/pdf/82049268.pdf

    // The steps of this method are as follows
    // 1. select a matrix (P) such that PAP^T has the smallest amount of "fill" elements when
    //    it is cholesky factorised. This is a NP-complete problem and hence non-trivial.
    //    There are many heuristic algorithms that find a P such that the fill in the cholesky is low
    //    e.g. nested disection algorithm & minimum degree algorithm
    //    For this implementation this step has not been completed (due to lack of time)
    //    This will mean that the (computationally expensive) symbolic factorisation and numeric
    //    factorisation is doing more iterations than is needed to still get a solution but it won't affect
    //    the accuracy of the result

    // 2. Symbolic Factorisation to determine the structure of the Cholesky lower diagonal matrix
    //    See the comments below for how this was implemented
    // 3. Numeric factorisation using the structure generated in the symbolic factorisation to calculate
    //    the values of the cholesky matrix
    // 4. Forward and Back substitution to solve the equation for x

    // makse sure the output is set to 0
    std::fill(this->output_sparse->begin(), this->output_sparse->end(), 0);

    // we need to create copies of the original matrix as it will be overwritten in place in the function
    // We also need a row and column vector to store the cholesky structure
    std::vector<T> values_original;
    std::vector<int> col_index_original;
    std::vector<int> col_index_cholesky;
    // as these are both the same length we can do the push back at the same time
    for (int i=0; i<nnzs; i++)
    {
        values_original.push_back(this->values->at(i));
        col_index_original.push_back(this->col_index->at(i));
        col_index_cholesky.push_back(this->col_index->at(i));
    }

    std::vector<int> row_position_original;
    std::vector<int> row_position_cholesky;
    for (int i=0; i<rows +1; i++)
    {
        row_position_original.push_back(this->row_position->at(i));
        row_position_cholesky.push_back(this->row_position->at(i));
    }
    
    // SYMBOLIC FACTORISATION
    // Start off with the two vectors: col_index_cholesky & row_position_cholesky
    // iterate through the col_index_cholesky and look backwards to the left to calculate the fill values
    // The structure of a column (j) in the cholesky matrix is given by the structure in the original matrix
    // at column (j) (hence why we take a copy to begin with)
    // This is combined with the structure of the columns in the cholesky matrix (L) whose 
    // off diagonal elements in row (j) are non-zero 
    // this involves in a given row (j) identifying if there are elements in the matrix that occur before the diagonal
    // the structure of these columns should be combined with the structure of the original matrix column (j)
    for (int r = 0; r < rows; r++)
    {
        int start_pair = row_position_cholesky[r]; // the starting CSR index of the elements in the current row
        int end_pair = row_position_cholesky[r+1]; // the end CSR index of the elements in the current row
        // if the 1st element in row is the diagonal then continue as there is no additional fill elements
        if (col_index_cholesky[start_pair] != r)  
        {
            // iterate over the columns that exist in the current row
            for (int j=start_pair; j< end_pair; j++)
            {
                int curr_col_idx = col_index_cholesky[j];
                // we only want columns that are less than the position of the diagonal (r) (looking backwards)
                if (curr_col_idx < r) 
                {
                    // curr_col_idx is the col index of one of the off diagonals in row (r)
                    // we want the structure of this column to be added to the structure of the current column (j)
                    int start_pair2 = row_position_cholesky[curr_col_idx];
                    int end_pair2 = row_position_cholesky[curr_col_idx +1];
                    // iterate over the elements in this column (making use of the fact that CSR and CSC are equal for Symmetric matrices)
                    for (int x=start_pair2; x<end_pair2; x++)
                    {
                        // we only want elements that are below the current row (as those that are above will be in the 
                        // upper part of the matrix when added to the current column (j)) (we will be deleting all upper
                        //  element later so this saves work later)
                        if (col_index_cholesky[x] > r)
                        {
                            // we want to insert the element col_index_cholesky[x] into the the current row in the col_index vector
                            // but first check we are not duplicating an existing element
                            // This may introduce more work than necessary as we need to remove duplicates anyway later
                            // so Possibly remove this condition to improve speed (but it result in more memory needed to store the duplicate elements)
                            if (!count(col_index_cholesky.begin() + row_position_cholesky[r], col_index_cholesky.begin() + row_position_cholesky[r+1], col_index_cholesky[x]))
                            {
                                auto newvec = col_index_cholesky.insert(col_index_cholesky.begin() + row_position_cholesky[r+1], col_index_cholesky[x]);
                                // as we are adding an element to the col_index vector we need to increment the row_position vector for all the subsequent elements
                                for (int z=r+1; z < rows+1; z++)
                                {
                                    row_position_cholesky[z] += 1;
                                }
                            }
                            // we want to insert the element (r) into the the row col_index_cholesky[x] in the col_index vector
                            // but first check we are not duplicating an existing element
                            // This may introduce more work than necessary as we need to remove duplicates anyway later
                            // so Possibly remove this condition to improve speed (but it result in more memory needed to store the duplicate elements)
                            if (!count(col_index_cholesky.begin() + row_position_cholesky[col_index_cholesky[x]], col_index_cholesky.begin() + row_position_cholesky[col_index_cholesky[x] + 1], r))
                            {
                                auto newvec2 = col_index_cholesky.insert(col_index_cholesky.begin() + row_position_cholesky[col_index_cholesky[x] + 1], r);
                                // as we are adding an element to the col_index vector we need to increment the row_position vector for all the subsequent elements
                                for (int z=col_index_cholesky[x]+1; z < rows+1; z++)
                                {
                                    row_position_cholesky[z] += 1;
                                }
                            }


                        }
                    }
                }
            }
        }
    }
    // remove the upper diagonal elements from the cholesky structure
    for (int r = 0; r < rows; r++)
    {
        int start_pair = row_position_cholesky[r];
        int end_pair = row_position_cholesky[r+1];
        int elem_removed = 0;
        for (int j=start_pair; j< end_pair; j++)
        {
            if (col_index_cholesky[j] > r)
            {
                col_index_cholesky.erase(col_index_cholesky.begin() + j);
                // this part of the code just allows us to avoid having the "for" loop to decrement 
                // the row_position as a nested loop and instead only update it once per cycle
                elem_removed += 1;
                end_pair -= 1;
                j-=1;
            }
        }
        // sort the column indices to be row major order and update the row_position values
        sort (col_index_cholesky.begin() + start_pair, col_index_cholesky.begin() + end_pair);
        for (int x = r+1; x< rows+1; x++)
        {
            row_position_cholesky[x] -= elem_removed;
        }
    }

    // remove the upper diagonal elements from original A matrix
    for (int r = 0; r < rows; r++)
    {
        int start_pair = row_position_original[r];
        int end_pair = row_position_original[r+1];
        int elem_removed = 0;
        for (int j=start_pair; j< end_pair; j++)
        {
            if (col_index_original[j] > r)
            {
                col_index_original.erase(col_index_original.begin() + j);
                values_original.erase(values_original.begin() + j);
                // this part of the code just allows us to avoid having the "for" loop to decrement 
                // the row_position as a nested loop and instead only update it once per cycle
                elem_removed += 1;
                end_pair -= 1;
                j-=1;
            }
        }
        // update the row position values according to how many elements were removed
        for (int x = r+1; x< rows+1; x++)
        {
            row_position_original[x] -= elem_removed;
        }
    }

    // Now we have created the col_index and row_position vectors for the cholesky Lower matrix structure
    // Now we want to take this structure and map the values from the original matrix onto this 
    // inserting "0" for the filled values
    for (int i =0; i< col_index_cholesky.size(); i++)
    {
        if (col_index_original[i] != col_index_cholesky[i])
        {
            col_index_original.insert(col_index_original.begin() + i, col_index_cholesky[i]);
            values_original.insert(values_original.begin() + i, 0);
        }
    }

    // NUMERIC FACTORISATION
    // For this we create a vector to store the values of the Lower Cholesky Matrix

    // Note the current state of our vectors are:
    // values_original = the original matrix values but with extra "0s" for the fill positions
    // row_position_cholesky = the correct row position indexing for the cholesky structure
    // col_index_cholesky = the correct col indices for the cholesky structure
    // Lower_cholesky = the vector to store the values of the numeric factorisation
    std::vector<T> Lower_Cholesky(values_original.size(), 0);
    
    // iterate over all the values in the filled sparse format (Lower Cholesky)
    // This is the same concept as the cholesky Dense fwd and back substitution but
    // now we are only iterating over the non-zero elements through the CSR format
    for (int i=0; i<rows; i ++)
    {
        int start_pair = row_position_cholesky[i];
        int end_pair = row_position_cholesky[i+1];

        // only iterating over the non-zero elements
        for (int j=start_pair; j<end_pair; j++)
        {
            T s = 0;
            if (col_index_cholesky[j] == i) // i.e we are on a diagonal
            {
                // use the formula for the diagonal element
                // sum of the square of the elements to the left of the diagonal in that row
                for (int k=start_pair; k<j; k++)
                {
                    s+= pow(Lower_Cholesky[k], 2);
                }
                // update the value of the L(j) element in the Cholesky decomposition matrix
                // L(j) = square root of the same value in the A mmatrix - the sum calculated above 
                Lower_Cholesky[j] = sqrt(values_original[j] - s);
            }

            else // we are not on the diagonal element 
            {
                // use the non diagonal element formula
                for (int k=start_pair; k<j; k++)
                {
                    //go given by the current column (col_idx[k])
                    for (int z = row_position_cholesky[col_index_cholesky[j]]; z< row_position_cholesky[col_index_cholesky[j] + 1]; z++)
                    {
                        // check if there is col_index_cholesky with the same value
                        if (col_index_cholesky[z] == col_index_cholesky[k])
                        {
                            // then take the product of the values and add it to s
                            s += Lower_Cholesky[z] * Lower_Cholesky[k];
                        }
                    }
                }
                Lower_Cholesky[j] = (values_original[j] - s) / Lower_Cholesky[row_position_cholesky[col_index_cholesky[j] + 1] - 1];
            }

        }
    }

    // FORWARD and BACK SUBSTITUTION
    // using Lower_cholesky, col_index_cholesky, row_position_cholesky, b

    // create a temporary vector y to store the values of the initial output
    std::vector<T> y (rows, 0);
    T s = 0;

    // iterate over the elements in the lower cholesky matrix and solve the
    // equation Ly=b and store the answer in y
    for (int i=0; i < rows; i++)
    {
        s=0;
        int start_pair = row_position_cholesky[i];
        int end_pair = row_position_cholesky[i+1];

        for (int k = start_pair; k < end_pair; k++)
        {
            s+=Lower_Cholesky[k] * y[col_index_cholesky[k]];
        }
        y[i] = (this->b->at(i) - s) / Lower_Cholesky[row_position_cholesky[i+1] -1];
    }
    this->output_sparse->at(rows - 1) = y[rows - 1] / Lower_Cholesky[row_position_cholesky[rows] - 1];

    // backward substitution (L^T)x = y (solving for x)
    // iterate over the elements in the transpose of the cholesky matrix (ie the upper values)
    for (int i=rows-2; i >=0; i--)
    {
        s=0;
        for (int ii = rows-2; ii >= i; ii--)
        {
            for (int jj = row_position_cholesky[ii + 2] - 1; jj >= row_position_cholesky[ii + 1]; jj--)
            {
                if (col_index_cholesky[jj] == i)
                {
                    s+=Lower_Cholesky[jj] * this->output_sparse->at(ii + 1);
                }
            }
        }
        this->output_sparse->at(i) = (y[i] - s) / Lower_Cholesky[row_position_cholesky[i+1] - 1];
    }

}

template <class T>
void CSRMat<T>::LU_Banded()
{
    // Here we are assuming a Banded Diagonal matrix
    // This has the implication that there is no additional values in the decomposed matrices 
    // i.e. the L & U structure (when combined) is identical to the original A matrix
    // This means we can avoid the symbolic factorisation stage
    // Unlike the cholesky we can solve using this method for Non-Symmetric matrices
    // makse sure the output is set to 0
    std::fill(this->output_sparse->begin(), this->output_sparse->end(), 0);
    int n = this->rows;
    int Lcounter = 0;
    int Ucounter = 0;
    int index_counter = 0;
    T s = 0;
    std::vector<T> y (n, 0); // to store temporary output before assigning it to the final result
    std::vector<T> Lrow (n+1, 0);
    std::vector<T> Urow (n+1, 0);

    // extract the structure of the Upper and Lower matrices from the CSR format
    for (int i=0; i< this->row_position->size() -1; i++)
    {
        int start_pair = this->row_position->at(i);
        int end_pair = this->row_position->at(i+1);
        for (int j=start_pair; j < end_pair; j++)
        {
            // test if we are on a diagonal
            if (this->col_index->at(j) == i)
            {
                Lcounter += 1;
                Ucounter +=1;
            }
            if (this->col_index->at(j) > 1)
            {
                Ucounter += 1;
            }
            else
            {
                Lcounter +=1;
            }
        }
        Lrow[i + 1] = Lcounter;
        Urow[i + 1] = Ucounter;
    }
    std::vector<T> Lvalues (Lcounter, 0);
    std::vector<T> Lcols (Lcounter, 0);
    std::vector<T> Uvalues (Ucounter, 0);
    std::vector<T> Ucols (Ucounter, 0);

    // Decompose the matrix into the Upper and Lower form using the structure generated above
    // iterate over the rows
    for (int i = 0; i < this->row_position->size() - 1; i++)
    {
        Lcounter = 0;
        Ucounter = 0;

        for (int j = 0; j < this->row_position->at(i+1) - this->row_position->at(i); j++)
        {
            // test if we are above a diagonal
            if (i > this->col_index->at(index_counter))
            {
                s = 0;
                for (int ii = 0; ii < this->col_index->at(index_counter); ii++)
                {
                    for (int jj = 0; jj < Urow[ii + 1] - Urow[ii]; jj++)
                    {
                        for (int kk = 0; kk < Lcounter; kk++)
                        {
                            // compute the actual value for the Lower matrix
                            if (Ucols[Urow[ii] + jj] == this->col_index->at(index_counter) 
                                && Lcols[Lrow[i + kk] == ii])
                            {
                                s += Lvalues[Lrow[i + kk]] * Uvalues[Urow[ii] + jj];
                            }
                        }
                    }
                }
                // Update the value in the stored Lower MAtrix
                Lvalues[Lrow[i] + Lcounter] = (this->values->at(index_counter) - s)/Uvalues[Urow[this->col_index->at(index_counter)]];
                Lcols[Lrow[i] + Lcounter] = this->col_index->at(index_counter);
                Lcounter += 1;
            }
            else
            {
                // for the diagonal element
                if (i == this->col_index->at(index_counter))
                {
                    Lvalues[Lrow[i] + Lcounter] = 1;
                    Lcols[Lrow[i] + Lcounter] = this->col_index->at(index_counter);
                }

                // For non diagonal elements in the upper matrix
                s = 0;
                for (int ii = 0; ii < i; ii++)
                {
                    for (int jj = 0; jj < Urow[ii + 1] - Urow[ii]; jj++)
                    {
                        for (int kk = 0; kk < Lcounter; kk++)
                        {
                            // compute the actual value for the upper matrix using the formula
                            if (Ucols[Urow[ii] + jj] == this->col_index->at(index_counter)
                            && Lcols[Lrow[i] + kk] == ii)
                            {
                                s += Lvalues[Lrow[i] + kk] * Uvalues[Urow[ii] + jj];
                            }
                        }
                    }
                }
                // Update the value in the stored Upper MAtrix
                Uvalues[Urow[i] + Ucounter] = this->values->at(index_counter) - s;
                Ucols[Urow[i] + Ucounter] = this->col_index->at(index_counter);
                Ucounter += 1;
            }
            
            index_counter += 1;
        }

        Lrow[i + 1] = Lcounter + Lrow[i] + 1;
        Urow[i + 1] = Ucounter + Urow[i];
    }

    // Forward substitution
    for (int i = 0; i < n; i++)
    {
        s = 0;
        for (int k = 0; k < Lrow[i + 1] - Lrow[i] - 1; k++)
        {
            s += Lvalues[Lrow[i] + k] * y[Lcols[Lrow[i] + k]];
        }
        // store values in y to use in backwwards substitution
        y[i] = this->b->at(i) - s;
    }

    // Backward substitution
    for (int i = n - 1; i >= 0; i--)
    {
        s = 0;
        for (int k = 0; k<Urow[i + 1] - Urow[i]; k++)
        {
            s += Uvalues[Urow[i] + k] * this->output_sparse->at(Ucols[Urow[i] + k]);
        }
        // store values in the output vector
        this->output_sparse->at(i) = (y[i] - s)/Uvalues[Urow[i]];
    }
}

// Matrix vector multiplication
// can be replaced by level 2 BLAS (such as SGEMV/DGEMV (general matrix-vector multiply) and
// SSYMV/DSYMV (if the matrix is symmetric))
// This is our implementation of this method
// we take advantage of the row ordering process to access elements that are near to each other 
// in contiguous memory in the inside loop
template <class T>
void CSRMat<T>::matA_vec_mult(std::vector<T> *vector_ptr, std::vector<T> *return_vec_ptr)
{
    // note the vector input and the resulting vector will both have length of "rows"
    // set the output vector to zero
    for (int i=0; i<this->rows; i++)
    {
        return_vec_ptr->at(i) = 0.0;
    }

    // loop over the rows
    for (int r=0; r<this->rows; r++)
    {
        // loop over the row pairing
        // std::cerr << "Print Row: " << r << "\t Values: ";
        for (int i=this->row_position->at(r); i < this->row_position->at(r+1); i++)
        {
            // std::cerr  << this->values[i] << " ," << this->col_index[i]  << "i=" << i << " r=" << r << "\t" ;
            return_vec_ptr->at(r) += this->values->at(i) * vector_ptr->at(this->col_index->at(i));
            
        }
        // std::cerr << std::endl;
    }
}

// Vector transpose - vector multiplication (dot product)
// can be replaced by level 1 BLAS SDOT, DDOT
template <class T>
T CSRMat<T>::vec_Trans_vec_mult(std::vector<T> *vector_left_ptr, std::vector<T> *vector_right_ptr)
{
    T return_val = 0;
    for (int i=0; i<this->rows; i++)
    {
        return_val += vector_left_ptr->at(i) * vector_right_ptr->at(i);
    }
    return return_val;
}

// This was never fully implemented an is hence commented out
// template <class T>
// void CSRMat<T>::Populate(T min_val, T max_val, int psd_factor, bool symmetric)
// {
//     // srand(time(0));
//     // this->symmetric = symmetric;
//     // this->max_val = max_val;
//     // this->min_val = min_val;
//     // this->psd_factor = psd_factor;

//     // for (int i=0; i<nnzs; i++)
//     // {
//     //     this->values[i] = min_val + rand() / (RAND_MAX / (max_val - min_val));
//     // }
//     // // allocate the row_position
//     // this->row_position[0] = 0; // first value is always 0
//     // for (int i=1; i<this->rows; i++)
//     // {
//     //     // fill rest of values with random numbers between 0 and nnzs not including nnzs
//     //     // these have to be in incrementing order hance we adda random amount to the previous
//     //     // with a max bound set by the nnzs
//     //     this->row_position[i] = std::min(this->row_position[i-1] + rand() % this->cols, this->nnzs - 1);
//     // }
//     // this->row_position[this->rows] = this->nnzs; // last value is always nnzs
    

//     // //allocate the col_index

//     // // create an array that contains all the possible column indices
//     // auto *possible_col_values = new int[this->cols] {0};
//     // for (int i=0; i<this->cols; i++)
//     // {
//     //     possible_col_values[i] = i;
//     // }
    
//     // for (int r=1; r<this->rows +1 ; r++)
//     // {
//     //     // iterate over the values between each pairining in the row_position
//     //     for (int i=this->row_position[r-1]; i<this->row_position[r]; i++)
//     //     {
//     //         // store the value of the number of non-zero values in the current row
//     //         int no_values_in_row = this->row_position[r] - this->row_position[r-1];
//     //         auto *temp_array = new int[no_values_in_row];
//     //         // create a random set of unique numbers from the potential column indices
//     //         // of length = the number of non-zero elements in the row
//     //         // assign these col_indices to the correct sace in the col_index array
//     //         std::sample(possible_col_values, possible_col_values + this->cols, temp_array,
//     //                 no_values_in_row, std::mt19937{std::random_device{}()});

//     //         for (int z=0; z< no_values_in_row; z++)
//     //         {
//     //             this->col_index[r-1 + z] = temp_array[z];
//     //         }
//     //         delete[] temp_array;
            
//     //     }
//     // }
//     // // make the first element 0
//     // // this->col_index[0] = 0; // can't do this as it would introduce a duplicate
//     // delete[] possible_col_values;
// }