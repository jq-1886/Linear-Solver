#include <iostream>
#include "DMat.h"
#include <memory>
#include <ctime>
#include <math.h>
#include <cstdlib>
#include <cmath>

template <class T>
DMat<T>::DMat(int rows, int cols, bool usr_provided_input_output)
{
    // set the values of this instance of the class as the values defined in the header
    this->rows = rows;
    this->cols = cols;
    this->no_values = rows*cols;
    this->usr_provided_input_output = usr_provided_input_output;

    // Check of the user provided no input and if so instantiate the values, b & output arrays
    // with the correct size
    if (!this->usr_provided_input_output)
    {
        // create the correct space for the values pointer (previously it was a null ptr)
        this->values.reset(new T[rows*cols]);
        // create the correct space for the b pointer (previously it was a null ptr)
        this->b.reset(new T[rows]);        
        // create the correct for the output of the x vector (previously was null ptr)
        this->output_dense.reset(new T[rows]);
    }
}

template <class T>
DMat<T>::DMat(int rows, int cols, std::shared_ptr<T[]> values_ptr, std::shared_ptr<T[]> b_ptr, std::shared_ptr<T[]> output_dense_ptr)
{
    // This is the overloaded constructor for this class and allows the user to directly instantiate
    // the class with their own values using shared pointers
    this->rows = rows;
    this->cols = cols;
    this->no_values = rows*cols;
    this->usr_provided_input_output = usr_provided_input_output;

    // set the values for values, b & output to the pointers given by the user
    this->values = values_ptr;
    this->b = b_ptr;
    this->output_dense = output_dense_ptr;

}

template <class T>
DMat<T>::~DMat()
{
    // as we are using shared pointers once they go out of scope all items stored / created should be deleted
    // hence the destructor is empty
}

template <class T>
void DMat<T>::Populate(T min_val, T max_val, int psd_factor, bool symmetric, int seed)
{
    // initialise the random seed using a user parameter
    srand(seed);
    this->symmetric = symmetric;
    this->max_val = max_val;
    this->min_val = min_val;
    this->psd_factor = psd_factor;

    // generate a dense matrix with random values following row major ordering
    for (int i=0; i < this->rows; i++)
    {
        for (int j=0; j< this->cols; j++)
        {
            this->values[i*this->rows + j] = min_val + rand() / (RAND_MAX / (max_val - min_val));
        }
    }

    // To create a symmetric matrix just adapt the one created above
    if (this->symmetric)
    {
        // iterate over the upper triangle and set values equal to lower triangle values
        int column_counter = 1;
        for (int i=0; i < this->rows; i++)
        {
            for (int j=column_counter; j<this->cols; j++)
            {
                this->values[i*this->rows + j] = this->values[j*this->rows + i];
            }
            column_counter += 1;
        }          
    }

    // Iterate over the diagonal and make positive definite by increasing 
    // the value and taking the absolute value
    for (int i=0; i < this->rows; i++)
    {
        this->values[i*this->rows +i] = abs((min_val + rand() / (RAND_MAX / (max_val - min_val))) * psd_factor);
    }
    
    
    // Assign random values to the b vector
    for (int i=0; i < this->rows; i++)
    {
        this->b[i] = min_val + rand() / (RAND_MAX / (max_val - min_val));
    }
}

template <class T>
void DMat<T>::PrintMatrix_A()
{
    std::cerr << "Printing Matrix A: " << std::endl;
    for (int i=0; i < this->rows; i++)
    {
        for (int j=0; j<this->cols; j++)
        {
            std::cerr << this->values[i *this->cols + j] << ", ";
        }
        std::cerr << std::endl;
    }
    std::cerr << std::endl;
}

template <class T>
void DMat<T>::PrintVector_b()
{
    std::cerr << "Printing Vector b: " << std::endl;
    for (int i=0; i < this->rows; i++)
    {
        std::cerr << this->b[i] << ", ";
    }
    std::cerr << std::endl << std::endl;
}

template <class T>
void DMat<T>::PrintVector_x()
{
    std::cerr << "Printing Vector x: " << std::endl;
    for (int i=0; i < this->rows; i++)
    {
        std::cerr << this->output_dense[i] << ", ";
    }
    std::cerr << std::endl << std::endl;
}

template <class T>
void DMat<T>::jacobi(double tolerance)
{
    // reset the vector x (this->output_dense) to 0
    this->output_dense.reset(new T[rows] {0});

    // check matrix A is square, if not end
    if (this->rows != this->cols)
    {
        std::cerr << "Matrix A is not square, pass in a square matrix" << std::endl;
        return;
    }

    T *vector_x_old = nullptr;
    vector_x_old = new T[this->rows] {0.0};
    // create a sum_holder variable to store the sum of the nondiagonal values
    T sum_holder = 0.0;   

    // create a variable to store the number of iterations until tolerance is reached
    int iterations = 0;
    // create a variable to flag to the do loop that the tolerance conditions have been satisfied
    int flag = 0;

    // loop to implement jacobi's method
    do
    {
        // move in row-major order
        for(int i = 0; i < this->rows; i++)
        {
            // revert sum_holder to 0 at the start of each loop
            sum_holder = 0;
            for(int j = 0; j < this->cols; j++)
            {
                // skip the sum_holder operation if on the main diagonal
                if(i == j) continue;
                {
                    // perform the summation part of Jacobi's method using the vector_x_old
                    // element a_{i,j} * x_{j}
                    sum_holder += this->values[j + i * this->cols] * vector_x_old[j];
                }
            }
             // finish calculating that step of the jacobi which is the b vector ith element - the sum_holder
             // all multiplied by 1 / diagonal element in row (i)
            this->output_dense[i] = (1 / this->values[i + i * this->cols]) * (this->b[i] - sum_holder);
            // compare with previous value and check if below the tolerance
            if (abs(vector_x_old[i] - this->output_dense[i]) < tolerance)
            {
                // if tolerance for an element in the x vector has been reached, increment the flag
                flag++;
            }
            // copy values from the newly populated x vector (this->output_dense) to the x_old to use 
            // them in the next iteration
            vector_x_old[i] = this->output_dense[i];
        }
        iterations ++;
        // check all elements have reached the desired tolerance, if not continue
    }while(flag < this->rows);
    // print to console using standard error channel the number of iterations performed
    std::cerr << "Iterations: " << iterations;
    // Free the storage on the heap
    delete[] vector_x_old;
}

template <class T>
void DMat<T>::gauss_seidel(double tolerance)
{
    // reset the vector x (this->output_dense) to 0
    this->output_dense.reset(new T[rows] {0});
    // check matrix A is square, if not end
    if (this->rows != this->cols)
    {
        std::cerr << "Matrix A is not square, pass in a square matrix" << std::endl;
        return;
    }
    // create a sum_holder variable to store the sum of the nondiagonal values
    T sum_holder = 0.0;
    // create a variable to store the number of iterations until tolerance is reached 
    int iterations =0;
    // create a variable to flag to exit the loop when the tolerance conditions have been satisfied
    int flag = 0;
    // create a variable to hold the ith element of the previous iteration
    T copy_previous_guess = 0;

    // loop to implement Guass-Seidel method
    while (flag < this->rows)
    {
        // move in row-major order
        for (int i=0; i<this->rows; i++)
        {
            // revert sum_holder to 0 at the start of each loop
            sum_holder = 0.0;
            // copy across the old values of vector x to compare for flagging purposes
            copy_previous_guess = this->output_dense[i];
            
            for (int j=0; j<this->cols; j++)
            {
                if (j == i) continue;
                {
                    // perform the summation part of gauss-seidel method using 
                    // (this->output_dense) as the previous x vector value
                    sum_holder += this->values[j + i * this->cols] * this->output_dense[j];
                }
            }
            // update the output x vector with the gauss_seidel formula
            this->output_dense[i] = (1 / this->values[i + i * this->cols]) * (this->b[i] - sum_holder);
            // compare with previous value and check if below the tolerance
            if (abs(copy_previous_guess - this->output_dense[i]) < tolerance)
            {
                // if tolerance for an element in the x vector has been reached, increment the flag
                flag++;
            }
        }
        iterations ++;
    }
    std::cerr << "Iterations: " << iterations;
}

template <class T>
void DMat<T>::successive_over_relaxation(double tolerance, double r_factor, double user_initial_guess)
{
    // reset the vector x (this->output_dense) to 0
    this->output_dense.reset(new T[this->rows] {0});
    // check matrix A is square
    if (this->rows != this->cols)
    {
        std::cerr << "Matrix A is not square, pass in a square matrix" << std::endl;
        return;
    }
    // create an array to hold the values of the x vector
    T *matrix_x = nullptr;
    matrix_x = new T[this->rows] {user_initial_guess};

    T sum_holder = 0.0;   
    int iterations = 0;
    // create a variable to flag to exit the loop when the tolerance conditions have been satisfied
    int flag = 0;
    // loop to implement successive over relaxation method
    do
    {
        // move in row-major order
        for(int i = 0; i < this->rows; i++)
        {
            // revert sum_holder to 0 at the start of each loop
            sum_holder = 0;

            for(int j = 0; j < this->cols; j++)
            {
                // if condition to skip the sum_holder operation if on the main diagonal
                if(i == j) continue;
                {
                    // perform the summation part of successive over relaxation method 
                    // using the previous guesses x_vector values
                    sum_holder += this->values[j + i * this->cols] * matrix_x[j];
                }
            }
            // perform the population of the new x vector (this->output_dense) according to the SOR formula
            this->output_dense[i] = this->output_dense[i] + r_factor * ((this->b[i] - sum_holder) / this->values[i + i * this->cols] - this->output_dense[i]);
            // compare with previous value and check if below the tolerance
            if (abs(matrix_x[i] - this->output_dense[i]) < tolerance)
                flag++;
            matrix_x[i] = this->output_dense[i];
        }
        iterations ++;
    }while(flag < this->rows);
    std::cerr << "Iterations: " << iterations; // Print iterations
    delete[] matrix_x; // Free the storage on the heap
}

template <class T>
void DMat<T>::LU_decomp()
{
    // This solver uses row major ordering and puts all elements of the lower
    // upper matrix into one matrix LU.

    // reset the vector x (this->output_dense) to 0
    this->output_dense.reset(new T[rows] {0});
    // check matrix A is square
    if (this->rows != this->cols)
    {
        std::cerr << "Matrix A is not square, pass in a square matrix" << std::endl;
        return;
    }

    // Decomposing Matrix A into upper, lower matrices 
    // create a Templated array to store temporary values (will be deleted at end of method)
    auto *LU = new T[this->rows * this->cols];
    T s = 0;
    int n = this->rows;

    // Doolittle algorithm (take the lower matrix diagonal elements = 1)
    for (int k = 0; k < n; k++)
    {
        // Upper Triangle iteration
        for (int i = k; i < n; i++)
        {
            // initialise the sum
            s = 0;
            // Take the sum of LU(i, p) * LU(p, k) over all p < k
            for (int p = 0; p < k; p++)
            {
                s += LU[i*n + p] * LU[p*n + k];
            }
            // this gives us the value for position LU(i, k) = A(i,k) - sum
            LU[i*n + k] = this->values[i*n + k] - s;
        }

        // Lower Triangle iteration
        for (int j = k+1; j < n; j++)
        {
            // initialise the sum
            s = 0;
            //Take the sum of LU(k,p) * LU(p,j)
            for (int p= 0; p < k; p++)
            {
                s += LU[k*n + p] * LU[p*n + j]; ///////
            }
            if (LU[k*n + k] == 0)
            {
                // This if statement prevents division by zero
                // in the following step and stops the iteration
                std::cerr << "Matrix A cannot be decomposed into LU components" << std::endl;
                return;
            }
            //this gives us the value for position LU(k,j) = (A(k,j) - sum)/diagonal (LU(k,k))
            LU[k*n + j] = (this->values[k*n + j] - s) / LU[k*n +k];
        }
    }

    // Forward substitution Ly = b (solving for y)
    // can use BLAS/LAPACK: (1) to find inverse of L (LAPACK: dgetrf) and 
    //                      (2) perform matrix-vector multiplication L^(-1) b to solve for y
    //                          (using triangle matrix vector multiplication e.g.DTRMV as inverse of
    //    

    // create a Templated array to store temporary values (will be deleted at end of method)
    auto *y = new T[n];
    // iterate forwards over the LU array using the Lower elements only and solve 
    // the equation Ly=b and update y as we go
    for (int i=0; i<n; i++)
    {
        s = 0;
        for (int k=0; k<i; k++)
        {
            s += LU[i*n + k] *y[k];
        }
        y[i] = (this->b[i] - s) / LU[i*n + i];
    }

    // backward substitution Ux = y (solving for x)
    // similar implementation of BLAS/LAPACK described in forward substitution can be used here
    // iterate backwards over the LU array (from last row)using the Upper elements only and solve the equation
    // Ux = y (using the y calculated above)
    for (int i = n-1; i>=0; i--)
    {
        s = 0;
        for (int k = i+1; k < n; k++)
        {
            s += LU[i*n + k] * this->output_dense[k];
        }
        // update the output array with the final solution
        this->output_dense[i] = (y[i] - s);
    }
    // delete any storage created during method
    delete[] LU;
    delete[] y;
}


template <class T>
void DMat<T>::Cholesky()
{
    // reset the vector x (this->output_dense) to 0
    this->output_dense.reset(new T[rows] {0});
    // check matrix A is square
    if (this->rows != this->cols)
    {
        std::cerr << "Matrix A is not square, pass in a square matrix" << std::endl;
        return;
    }

    // Cholesky Decomposition A = LL* (assuming A is positive definite)
    // Decomposing Matrix A into upper, lower matrices 
    // create a Templated array to store temporary values (will be deleted at end of method)
    auto *L = new T[this->rows * this->cols];
    T s = 0;
    int n = this->rows;

    L[0] = sqrt(this->values[0]);
    // iterate over the Lower elements in the matrix A
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < i+1; j++)
        {
            s = 0;
            if (i == j) // checking if we are on the diagonal element
            {
                // use the formula for the diagonal element
                // sum of the square of the elements to the left of the diagonal in that row
                for (int k = 0; k < j; k++)
                {
                    s+= pow(L[j*n + k], 2);
                }
                if (this->values[j*n + j] - s < 0)
                {
                    // This prevents square rooting a negative number in the
                    // following step nd stops the iteration
                    std::cerr << "Matrix A cannot be Cholesky decomposed: Matrix not positive definite" << std::endl;
                    return;
                }
                // update the value of the L(j,j) element in the Cholesky decomposition matrix
                // L(j,j) = square root of the same value in the A mmatrix - the sum calculated above 
                L[j*n + j] = sqrt(this->values[j*n + j] - s);
            }
            else // we are not on the diagonal element 
            {
                if (L[j*n + j] == 0)
                {
                    // This prevents dividing by zero
                    std::cerr << "Matrix A cannot be Cholesky decomposed: Matrix not positive definite" << std::endl;
                    return;
                }
                // use the non diagonal element formula
                for (int k = 0; k < j; k++)
                {
                    // look to the left of the current element in row (i)
                    // take to the sum of all these elements multiplied by their corresponding 
                    // element in the same column but in row (j) 
                    s += L[i*n + k]*L[j*n + k];
                }
                // update the value of cholesky matrix
                L[i*n + j] = (this->values[i*n + j] - s) / L[j*n + j];
            }
        }
    }

    // forward substitution Ly = b (solving for y)
    // Forward substitution Ly = b (solving for y)
    // can use BLAS/LAPACK: (1) to find inverse of L (LAPACK: dgetrf) and 
    //                      (2) perform matrix-vector multiplication L^(-1) b to solve for y
    //                          (using triangle matrix vector multiplication e.g.DTRMV as inverse of
    //    
    auto *y = new T[n];
    // iterate over the elements in the lower cholesky matrix and solve the
    // equation Ly=b and store the answer in y
    for (int i = 0; i < n; i++)
    {
        s = 0;
        for (int k = 0; k < i; k++)
        {
            s += L[i*n + k] * y[k];
        }
        y[i] = (this->b[i] - s) / L[i*n + i];
    }

    // backward substitution (L^T)x = y (solving for x)
    // iterate over the elements in the transpose of the cholesky matrix (ie the upper values)
    for (int i = n-1; i >=0; i--)
    {
        s = 0;
        for (int k = i+1; k < n; k++)
        {
            s += L[k*n + i] * this->output_dense[k];
        }
        // update the output array with the final solution
        this->output_dense[i] = (y[i] - s) / L[i*n + i];
    }
    // delete storage generated during the method
    delete[] L;
    delete[] y;
}

template <class T>
void DMat<T>::Conjugate_Grad(double tolerance)
{
    // https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
    // algorithm follws that given on page 50 of the above

    // reset the vector x (this->output_dense) to 0
    this->output_dense.reset(new T[rows] {0});
    // check matrix A is square
    if (this->rows != this->cols)
    {
        std::cerr << "Matrix A is not square, pass in a square matrix" << std::endl;
        return;
    }
    // initialise Templated arrays and value storage for calculations
    auto *r = new T[this->rows] {0}; // residual vector
    auto *d = new T[this->rows] {0}; // the update to our guess
    T sigma_new = 0; // this is the value of the residual
    T sigma_0 = 0;
    auto *q = new T[this->rows] {0}; // this will be the product of A . d
    T alpha = 0; // this is the length of the step we will take
    T sigma_old = 0; // this is to store the original sigma at each loop
    T beta = 0; // this is the ratio of the new residual : old residual
    auto *temp2 = new T[this->rows] {0}; // this will be product of A . x

    // In this method we make repeated use of matrix vector multiplication and
    // vector dot product multiplication, this could be an area where using BLAS  / LaPAK calls
    // would improve efficiency. We have choosen in our implementation to write our own versions
    // of these functions to show the principal behind them instead of linking to a library.
    // However in a "real world" implementation we would link to the available libraries.

    // initialise "r" following r = b - Ax
    matA_vec_mult(this->output_dense.get(), r); // r = Ax
    for (int i=0; i<this->rows; i++)
    {
        r[i] = this->b[i] - r[i];
    }
    // initialise d by setting it to r
    for (int i=0; i<this->rows; i++)
    {
        d[i] = r[i];
    }

    // initialise sigma_new by setting it to r^T . r
    sigma_new = vec_Trans_vec_mult(r, r);

    // initialise sigma_0 with the value of sigma_new
    sigma_0 = sigma_new;

    
    int iterations = 0;
    while ((iterations < 500000) && (sigma_new > pow(tolerance, 2) * sigma_0))
    
    //while (iterations < 100)
    {
        // q = A . d
        matA_vec_mult(d, q);

        // calculate d^T . q and assign it to temp
        T temp = 0;
        temp = vec_Trans_vec_mult(d, q);
        // initialise alpha: alpha = sigma_new / temp
        alpha = sigma_new / temp;


        // update x with our new approximation
        // x = x + alpha * d
        for (int i=0; i<this->rows; i++)
        {
            this->output_dense[i] += alpha *d[i];
        }

        // update the residual 
        // this is using the more computationally expensive method
        // r = b - Ax
        // calculate Ax and assign it to temp2
        matA_vec_mult(this->output_dense.get(), temp2);
        for (int i=0; i<this->rows; i++)
        {
            r[i] = this->b[i] - temp2[i];
        }

        // update sigma_old with the current value of sigma
        sigma_old = sigma_new;

        // update sigma_new = r^T . r
        sigma_new = vec_Trans_vec_mult(r, r);

        // calculate beta = sigma_new / sigma_old
        beta = sigma_new / sigma_old;

        // update d = r + beta * d
        for (int i=0; i<this->rows; i++)
        {
            d[i] = r[i] + beta * d[i];
        }
        iterations ++;

    }
    std::cerr << "Iterations: " << iterations;

    // delete all new storage created on the heap during the method
    delete[] temp2;
    delete[] q;
    delete[] d;
    delete[] r;
}

// Matrix vector multiplication
// can be replaced by level 2 BLAS (such as SGEMV/DGEMV (general matrix-vector multiply) and
// SSYMV/DSYMV (if the matrix is symmetric))
// This is our implementation of this method
// we take advantage of the row ordering process to access elements that are near to each other 
// in contiguous memory in the inside loop
template <class T>
void DMat<T>::matA_vec_mult(T *vector_ptr, T *return_vec_ptr)
{
    for (int i=0; i<this->rows; i++)
    {
        T sum = 0;
        for (int j=0; j<this->cols; j++)
        {
            sum += this->values[i*this->cols + j] * vector_ptr[j];
        }
        return_vec_ptr[i] = sum;
    }
}

// Vector transpose - vector multiplication (dot product)
// can be replaced by level 1 BLAS SDOT, DDOT
template <class T>
T DMat<T>::vec_Trans_vec_mult(T *vector_left_ptr, T *vector_right_ptr)
{
    T return_val = 0;
    for (int i=0; i<this->rows; i++)
    {
        return_val += vector_left_ptr[i] * vector_right_ptr[i];
    }
    return return_val;
}


