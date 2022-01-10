# Run Instructions
To see all the solvers called run the following command from the main repository

g++ MainSolver.cpp DMat.h DMat.cpp CSRMat.h CSRMat.cpp -std=c++17
Then run the executeable file that is created

This will demonstrate all the solvers being called and solving sample matrices (including outputting whether they solved correctly)

# Using the Library
There are 2 classes included which provide a variety of solvers for Dense and Sparsely formatted Matrices. 
The following describes the structure of each of the classes and the methods avaailable to the user to solve
systems of linear equations of the form Ax = b

The classes are Templated so the user can decide which "Type" of variable they wish to instantiate the class with. Note the inputs must be of consistent type (e.g. the input matrix A the vector b and the output vector x should be of the same type)

## Dense Solvers: DMat

### User Instructions
Instantiate the class with user inputs e.g.:
~~~
auto dense_user_matrix = std::make_shared<DMat<double>>(rows, cols, user_values_ptr, user_RHS_ptr, user_output_ptr);
~~~  
#### DMat(int rows, int cols, std::shared_ptr<T[]> values_ptr, std::shared_ptr<T[]> b_ptr, std::shared_ptr<T[]> output_dense_ptr)  
Constructor for the Dense Matrix Class
*Parameters*  
rows = Integer number of rows in A matrix  
cols = Integer number of columns in A matrix  
values_ptr = Shared Pointer to a Templated array containing the matrix A stored in Row major ordering  
b_ptr = Shared Pointer to a Templated array containing the vector b  
output_dense_ptr = Shared Pointer to a Templated array of size (rows) to store the answer  

Alternatively the user can instantiate the class without any user input matrix e.g.:
~~~
auto dense_sym_matrix = std::make_shared<DMat<double>>(rows, cols, false);
~~~
#### DMat(int rows, int cols, bool usr_provided_input_output)  
Overloaded Constructor for the Dense Matrix Class
*Parameters*  
rows = Integer number of rows in A matrix  
cols = Integer number of columns in A matrix  
usr_provided_input_output = Bool (Set to false if the user does not wish to provide a System of equations to solve, Set to True if the user wishes just to instantiate the class and later set the values of the problem.) e.g. below:  
~~~
auto dense_user_matrix = std::make_shared<DMat<double>>(rows, cols, true)
dense_user_matrix->values = user_values_ptr
dense_user_matrix->b = user_RHS_ptr
dense_user_matrix->output_dense = user_output_ptr
~~~

### Methods available
#### void PrintMatrix_A()
Print the current A matrix by iterating over all the elements in matrix A and prints to the console row by row using the standard error channel.  

#### void PrintVector_b()
Print the current b vector by iterating over all the elements in vector b and prints to the console using the standard error channel.  

#### void PrintVector_x()
Print the current x output vector by iterating ver all the elements in vector x and prints to the console using the standard error channel.  

#### void Populate(T min_val, T max_val, int psd_factor, bool symmetric, int seed=0)  
If the user has specified ```user_provided_input_output = false``` when instantiating the class then the populate function can be called to generate a random sample problem of the size specified. The function first creates a random dense matrix and then ```if symmetric == true``` it will set the upper triangle values equal to the lower riangle values. Finally scales the diagonal values by the positive definite factor and returns the absolute value to ensure a positive definite matrix.  
*Parameters*  
min_val = Templated variable for the approx min val of the generated matrix  
max_val = Templated variable for the approx max val of the generated matrix  
psd_factor = Integer for the factor by which to increase the size of the elements on the main diagonal  
symmetric = Boolean to choose whether to generate a symmetric matrix or not  
seed = Integer to have the ability to change the seed of the random generator  

#### void jacobi(double tolerance=10e-10)
Solve current system using Jacobi Iterative Algorithm. Suitable for Symmetric Real Positive Definite Matrices. Not guaranteed to converge Stable for Non-Symmetric Matrices. This function has the following element-wise iterations  
<img src="https://render.githubusercontent.com/render/math?math=x_{i}^{k%2B1}%20=%20\frac{1}{\alpha_{ii}}\left(b_i%20-%20\sum_{j\neq%20i}\alpha_{ij}x_{j}^{k}\right),%20i%20=%201,2,...,n">  
where x^k is the kth approcimation of the vector x. It is noted that x^k and x^(k+1) cannot be overridden and so the minimum memory requirement of this method is of two 1xn vectors, where A is an nxn matrix.  
*Parameters*  
tolerance = Double for the user to set the value for the tolerance of the convergence that is acceptable  

#### void gauss_seidel(double tolerance=10e-10)
Solve current system using Gauss Seidel Iterative Algorithm. Suitable for Symmetric Real Positive Definite Matrices. Not guaranteed to converge Stable for Non-Symmetric Matrices. This function follows the Gauss-Seidel iterative method, which has the following element-wise iterations, that take advantage of forward substitution:  
<img src="https://render.githubusercontent.com/render/math?math=x_{i}^{k%2B1}%20=%20\frac{1}{\alpha_{ii}}\left(b_i%20-%20\sum_{j=i}^{i-1}\alpha_{ij}x_{j}^{k%2B1}%20-%20\sum_{j=i%2B1}^{n}\alpha_{ij}x_{j}^{k}\right),%20i%20=%201,2,...,n">  
This method is very similar to Jacobi’s and it is noted that the computation of xi^(k+1) uses the already computed values of xi^(k+1), and no holding array to store the previous iteration of vector x is needed  
*Parameters*  
tolerance = Double for the user to set the value for the tolerance of the convergence that is acceptable  

#### void successive_over_relaxation(double tolerance=10e-10, double r_factor=1.0)
Solve current system using Successive Over Relaxation Iterative Algorithm. Suitable for Symmetric Real Positive Definite Matrices. Not guaranteed to converge Stable for Non-Symmetric Matrices as this is a variant of the Gauss Seidel Method. This function follows the successive over relaxation “SOR” iterative method, which has the following element-wise iterations, that take advantage of forward substitution:  
<img src="https://render.githubusercontent.com/render/math?math=x_{i}^{k%2B1}%20=%20(1%20-%20\omega)x_{i}^{k}\frac{\omega}{\alpha_{ii}}\left(b_i%20-%20\sum_{j<i}\alpha_{ij}x_{j}^{k%2B1}%20-%20\sum_{j>i}\alpha_{ij}x_{j}^{k}\right),%20i%20=%201,2,...,n">  
The similarity between SOR and Gauss-Seidel is noted in that Gauss-Seidel is essentially SOR but with a fixed relaxation factor of 1. It is also noted that in the element-wise iterations, the population of the vector x can be re-written as  <img src="https://render.githubusercontent.com/render/math?math=x_{i}%20%2B%20\omega\left(\frac{b_{i}%20-%20\sigma}{\alpha_{ii}}%20-%20x_{i}\right),%20where%20\sigma"> is the sum product of the non-diagonal elements of matrix A and vector x, thus saving one multiplicative step in each iteration.   
*Parameters*  
tolerance = Double for the user to set the value for the tolerance of the convergence that is acceptable  
r_factor = Double Set by the user to determine the rate (If A is symmetric positive definite then 0 < r < 2)  

#### LU_decomp()
Solve current system using LU decomposition. Suitable for Symmetric & non-Symmetric Real Positive Definite Matrices. The function first factorizes the input matrix A = LU using Dolittle’s algorithm, and solves the equation Ax = b using forward and backward substitution, the result x is then stored in the output_dense object of the DMat class. Here, the function treats the lower and upper triangular matrices (L, U) as one matrix LU, with its upper triangular elements (including diagonal) equals U and lower triangular elements (not including diagonal) equals L. 
The function can return an error if A is not LU decomposable. While this method is slow and requires more storage due to the number of nested loops, it gives an explicit solution for x

#### Cholesky()
Solve current system using Cholesky decomposition. Suitable for Symmetric Real Positive Definite Matrices. The function first performs a Cholesky defactorisation on the input matrix ( A = LL*, where * denotes the transpose) and, similar to LU_decomp, it solves the equation Ax = b using forward and backward substitution then stores the result x in output_dense. The function calculates the lower triangular matrix L (defaulting to row-major-ordering) and treats L in column-major-ordering as L*.
This method gives an explicit solution for x

#### void Conjugate_Grad(double tolerance=10e-10)  
Solve current system using Conjugate Gradient method suitable for symmetric and non-symmetric 
real Positive Definite Matrices. Iterative method with a hardcoded upper limit of 500000 iterations.  
*Parameters*  
tolerance = Double for the user to set the value for the tolerance of the convergence that is acceptable  

#### void matA_vec_mult(T *vector_ptr, T *return_vec_ptr)
Compute the Matrix (A) vector multiplication for a dense matrix (This is instead of making a call to BLAS / LAPACK Level 2 which would be a quicker implementation but not done here to demonstrate how the algorithm works)   
*Parameters*  
vector_ptr = Pointer to Templated  array containing the vector for multiplication  
return_vec_ptr = Pointer to Templated  array for the output to be stored  

#### T vec_Trans_vec_mult(T *vector_left_ptr, T *vector_right_ptr);
Compute the dot product of 2 vectors. (This is instead of making a call to BLAS / LAPACK Level 1 which would be a quicker implementation but not done here to demonstrate how the algorithm works)  
*Parameters*  
vector_left_ptr = Pointer to Templated  array containing the left vector for multiplication  
vector_right_ptr = Pointer to Templated  array containing the vector for multiplication  
*Returns*
Templated variable with the dot product of the vectors  


## Sparse Solvers: CSRMat

### User Instructions
Instantiate the class with user inputs e.g.:  
~~~
auto *sparse_user_matrix = new CSRMat<long double>(rows1, cols1, nnzs1, user_values_ptr1, user_row_ptr1, user_col_ptr1, user_RHS_ptr1, user_output_ptr1);
~~~  
#### CSRMat(int rows, int cols, int nnzs, std::vector<T>* values_ptr, std::vector<int>* row_position_ptr, std::vector<int>* col_index_ptr, std::vector<T>* b_ptr, std::vector<T>* output_sparse_ptr)
Constructor for the CSR Matrix classs
*Parameters*  
rows = Integer number of rows in A matrix  
cols = Integer number of columns in A matrix  
nnzs = Integer number of number of non zero values in A matrix
values_ptr = Pointer to Templated vector (STL) that contains the non-zero values of the sparse matrix in CSR format  
row_position_ptr = Pointer to Templatedvector (STL) that contains the index in the values_ptr array where each new row starts and has the nnzs at the end  
col_index_ptr = Pointer to Templated vector (STL) that contains the column index in the matrix of each of the values given in the values_ptr  
b_ptr = Pointer to Templated vector (STL) that contains the vector b   
output_dense_ptr = Pointer to Templated vector (STL) of size (rows) to store the answer  

### Methods available
#### void PrintMatrix_A()
Print the current A matrix  in CSR format: outputs the 3 vectors for the values, row position & column index  

#### void PrintVector_b()
Print the current b vector  

#### void PrintVector_x()
Print the current x output vector  

#### void sparse2dense(T *dense_mat_ptr)
*Parameters*  
dense_mat_ptr = Pointer to a Templated vector (STL) to store the values of the matrix A in dense format in Row Major ordering

#### void jacobi(double tolerance=10e-10)
Solve current system using Jacobi Iterative Algorithm on CSR format. Suitable for Symmetric & non-Symmetric Real Positive Definite Matrices  
*Parameters*  
tolerance = Double for the user to set the value for the tolerance of the convergence that is acceptable  

#### void successive_over_relaxation(double tolerance=10e-10, double r_factor=1.0)
Solve current system using Successive Over Relaxation Iterative Algorithm on CSR format. Suitable for Symmetric Real Positive Definite Matrices. Not guaranteed to converge Stable for Non-Symmetric Matrices as this is a variant of the Gauss Seidel Method.  
*Parameters*  
tolerance = Double for the user to set the value for the tolerance of the convergence that is acceptable  
r_factor = Double Set by the user to determine the rate (If A is symmetric positive definite then 0 < r < 2)  

#### void gauss_seidel(double tolerance=10e-10)
Solve current system using Gauss Seidel Iterative Algorithm CSR format. Suitable for Symmetric Real Positive Definite Matrices. Not guaranteed to converge Stable for Non-Symmetric Matrices.  
*Parameters*   
tolerance = Double for the user to set the value for the tolerance of the convergence that is acceptable  

#### void Conjugate_Grad(double tolerance=10e-10)
Solve current system using Conjugate Gradient method on CSR format. Suitable for symmetric and non-symmetric 
real Positive Definite Matrices. Iterative method with a hardcoded upper limit of 500000 iterations.  
*Parameters*   
tolerance = Double for the user to set the value for the tolerance of the convergence that is acceptable  

#### void Cholesky()
Solve current system using Cholesky decomposition on CSR format. Suitable for Symmetric Real Positive Definite Matrices. This is not a a convergent solver so it solves explicitly in a set number of iterations depending on the size of the matrix.  

#### void LU_Banded()
Solve current system using LU decomposition on CSR format. Suitable for Banded Non-Symmetric & symmetric Real Positive Definite Matrices. The function performs an LU defactorisation on A into lower and upper matrices (L, U respectively) in CSR format using Doolittle’s algorithm and solves the equation Ax = b using forward and backward substitution. The result is then stored in the output_dense object of the DMat class. 
The sparse LU decomposition here can be performed on banded symmetric and non-symmetric matrix A. The function yields an exact answer of x.

#### void matA_vec_mult(std::vector<T> *vector_ptr, std::vector<T> *return_vec_ptr)
Compute the Matrix (A) vector multiplication for a dense matrix (This is instead of making a call to BLAS / LAPACK Level 2 which would be a quicker implementation but not done here to demonstrate how the algorithm works)  
*Parameters*   
vector_ptr = Pointer to Templated vector (STL) containing the vector for multiplication  
return_vec_ptr = Pointer to Templated vector (STL) for the output to be stored   

#### T vec_Trans_vec_mult(T *vector_left_ptr, T *vector_right_ptr);
Compute the dot product of 2 vectors. (This is instead of making a call to BLAS / LAPACK Level 1 which would be a quicker implementation but not done here to demonstrate how the algorithm works)  
*Parameters*   
vector_left_ptr = Pointer to Templated vector (STL) containing the vector for multiplication    
vector_right_ptr = Pointer to Templated vector (STL) for the output to be stored   
*Returns*  
Templated variable with the dot product of the vectors  


# Example of code for user to run Solvers on user defined matrices:
## Dense Solver
~~~
the user values stored somewhere
int rows = 9;
int cols = 9;
double user_values[] = {20, 0, 0, 0, 1, 0, 4, 0, 0,
                        0, 20, 0, 0, 1, 0, 0, 0, 2,
                        0, 0, 20, 0, 0, 1, 8, 0, 0,
                        0, 0, 0, 20, 0, 6, 0, 0, 0,
                        1, 1, 0, 0, 20, 0, 0, 1, 0,
                        0, 0, 1, 6, 0, 20, 0, 5, 0,
                        4, 0, 8, 0, 0, 0, 20, 1, 0,
                        0, 0, 0, 0, 1, 5, 1, 20, 1,
                        0, 2, 0, 0, 0, 0, 0, 1, 20};
double user_RHS[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
~~~
Shared pointers to arrays (cannot use make_shared with arrays, enforce c++17 standard so that default deleter for shared pointer to arrays is present delete[])
~~~
std::shared_ptr<double[]> user_values_ptr(new double[rows*cols]); //need to try do this without "new" not possible because of pointer to array does not support "make_shared"
for (int i=0; i<rows*cols; i++)
{
    user_values_ptr[i] = user_values[i];
}

std::shared_ptr<double[]> user_RHS_ptr(new double[rows]); //need to try do this without "new"
for (int i=0; i<rows; i++)
{
    user_RHS_ptr[i] = user_RHS[i];
}
std::shared_ptr<double[]> user_output_ptr(new double[rows]);
~~~
Instantiate the Class with a shared pointer
~~~
auto dense_user_matrix = std::make_shared<DMat<double>>(rows, cols, user_values_ptr, user_RHS_ptr, user_output_ptr);
~~~
Call the functions using the class pointer
~~~
dense_user_matrix->jacobi()
~~~

