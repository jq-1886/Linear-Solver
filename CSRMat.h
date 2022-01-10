#pragma once
#include <memory>
#include <vector>

// This file contains the Sparse matrix class which instantiates a Sparse Matrix in CSR format
// and includes a number of linear solvers which are suited to a variety of problems
// The input matrix is assumed to be Real Positve Definite for all solvers
// In addition some solvers have further restrictions on the types of input matrix

// Matrices are stored in Compressed Row Format in vectors. We choose to use vectors in the 
// case of a sparse solver so as to avoid mistakes in handling memory as containers are a safer 
// approach to handling vectors where the lenght can change. These are accessed through pointers
// so as to reduce the need to copy the vectors and take up additional memory

// The following methods all take a pointer to the vecotors from the user as an input
// We also attempted to build a Sparse matrix generator but decided it was better to just
// provide sample matrices in the main file (simulating a user input).

// Some additional heap storage is allocated within methods but this is all destroyed before the method ends

// In choosing the structure for the class we decided to include all the solvers within the class
// where the matrix is constructed. It would be possible to create each solver in its own seperate
// sub-class that inheirits the structure and values of the matrix from the main Dense class
// However we decided that this was merely implementing hierarchy for the sake of it as all the solvers need access
// access to all the elements that are constructed in this Parent class (CSRMat). This is a good indication
// that the solvers should actually be in the same class as where the items they need are constructed.
// Due to the volume of solvers it could be argued that some solvers could be clearer if in seperate files but this
// is purely an asthetic consideration that we did not see the benefit of (admittedly for a commercial
// linear solver package that has far more solver, each in more detail than presented here this is a more 
// important consideration to enhance the sustainability of the code)

// There are occassional calculations that are shared between different solvers but again we did not see
// this as sufficiently large overlap to warrant creating subclasses however this could be implemented.


template <class T>
class CSRMat
{
public:
    // constructor
    CSRMat(int rows, int cols, int nnzs, bool usr_provided_input_output);

    // Overloaded constructor to allow the user to instantiate the Sparse Matrix class directly
    // The user must pass in pointer to a vector of Templated Type In CSR format:
    // Values vector (Templated)
    // Col index vector (Integers)
    // Row position vector(Integers)
    // b vector (Templated)
    // output vector (Templated)
    CSRMat(int rows, int cols, int nnzs, std::vector<T>* values_ptr, std::vector<int>* row_position_ptr, std::vector<int>* col_index_ptr, std::vector<T>* b_ptr, std::vector<T>* output_sparse_ptr);

    // destructor
    virtual ~CSRMat();

    // print matrix (A contained inside class)
    void PrintMatrix_A();

    // Print vector (b contained inside class)
    void PrintVector_b();

    // Print general
    void PrintVector_x();

    // sparse to dense converter
    // takes a pointer to an array of size rows*cols as an input and populates this with
    // the dense version of the "A" matrix currently stored in CSR form
    void sparse2dense(T *dense_mat_ptr);

    // populate with a randomly generated matrix
    // the code for this is incomplete and hence commented out. Instead of generating
    // matrices automatically we have included several test matrices in the maintest file
    // void Populate(T min_val, T max_val, int psd_factor, bool symmetric);

    // solvers
    // Jacobi: Symmetric & Non-Symmetric
    void jacobi(double tolerance=10e-10);

    // SOR: Symmetric (not guaranteed to converge for all Non-Symmetric Matrices)
    void successive_over_relaxation(double tolerance=10e-10, double r_factor=1.0);

    // G-S: Symmetric (not guaranteed to converge for all Non-Symmetric Matrices)
    void gauss_seidel(double tolerance=10e-10);

    // CJ: Symmetric & Non-Symmetric
    void Conjugate_Grad(double tolerance=10e-10);

    // Cholesky: Symmetric only implementation
    void Cholesky();

    // LU Decomp: Banded Symmetric Matrix
    void LU_Banded();

    // methods
    // matrix vector multiplication
    void matA_vec_mult(std::vector<T> *vector_ptr, std::vector<T> *return_vec_ptr);
    T vec_Trans_vec_mult(std::vector<T> *vector_left_ptr, std::vector<T> *vector_right_ptr);

    // These are the pointers to the different vectors -- shared pointers were not recommended for use with
    // vectors after doing research. This leaves more risk for memory to leak if not deleted from the heap
    // However we have explicitly required the user to provide these pointer and no other memory is created
    // on the heap so there should be no memory leak. The user must provide the pointer to the output
    // so they have ownership of that memory already and are aware of it!
    // values (non zero values)
    std::vector<T>* values = nullptr;

    // row position (length rows + 1)
    std::vector<int>* row_position = nullptr;

    // column index for each non zero value
    std::vector<int>* col_index = nullptr;

    // b (vector b RHS) (length rows of LHS)
    std::vector<T>* b = nullptr;

    // output_sparse (the returned x vector)
    std::vector<T>* output_sparse = nullptr;

    // rows
    int rows = -1;
    int cols = -1;
    int nnzs = -1;

    // initialise the min and max values for the random ssample matrix 
    // (can be changed by user later)
    T min_val = -10;
    T max_val = 10;
    int psd_factor = 10;

protected:
    // preallocated
    bool usr_provided_input_output = false;
    bool symmetric = true;

private:
    // size of values
    int no_values = -1;

};