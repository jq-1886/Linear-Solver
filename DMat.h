#pragma once
#include <memory>

// This file contains the Dense matrix class which instantiates a Dense Matrix in Row Major ordering
// and includes a number of linear solvers which are suited to a variety of problems
// The input matrix is assumed to be Real Positve Definite for all solvers
// In addition some solvers have further restrictions on the types of input matrix

// Matrices are stored in Row Major ordering using arrays,

// To reduce memory space required the arrays are all accessed through shared pointers which
// removes (most but not all!) of the need to worry about memory leak
// This is particularly improtant for returning matrices to users and this ensures that when the class
// is destroyed all heap storage of the matrices will be deleted

// More sophisticated templated containers that can be manipulated more easily
// such as vectors, are more suited to storing these matrices.
// However for the purposes of this assignment to demonstrate our understanding of how to implement
// fundamental C++ prinipals we choose to use arrays for the dense matrix solvers
// As dense matrix solvers never need to adjust the elements once inputtedd (no increasing / decreasing
// of the length) this makes arrays very suitable. When the ength of arrays needs to be adjusted it becomes
// more efficient to implement a container in the STL such as vectors as this will lead to better processing
// and less of a chance to accidently make errors in handling memory 

// To demonstrate our understanding of how to implement the STL we have built our Sparse matrix solver using
// vectors instead of arrays for the reasons outlined above.

// Some additional heap storage is allocated within methods but this is all destroyed before the method ends

// In choosing the structure for the class we decided to include all the solvers within the class
// where the matrix is constructed. It would be possible to create each solver in its own seperate
// sub-class that inheirits the structure and values of the matrix from the main Dense class
// However we decided that this was merely implementing hierarchy for the sake of it as all the solvers need access
// access to all the elements that are constructed in this Parent class (DMat). This is a good indication
// that the solvers should actually be in the same class as where the items they need are constructed.
// Due to the volume of solvers it could be argued that some solvers could be clearer if in seperate files but this
// is purely an asthetic consideration that we did not see the benefit of (admittedly for a commercial
// linear solver package that has far more solver, each in more detail than presented here this is a more 
// important consideration to enhance the sustainability of the code)

// There are occassional calculations that are shared between different solvers but again we did not see
// this as sufficiently large overlap to warrant creating subclasses however this could be implemented.

template <class T>
class DMat
{
public:
    // constructor
    DMat(int rows, int cols, bool usr_provided_input_output);

    // Overloaded constructor to allow the user to instantiate the Dense Matrix class directly
    // The user must pass in a shared pointer to an array of Templated Type for the:
    // Values Matrix (stored in Row major ordering)
    // b vector
    // output vector
    DMat(int rows, int cols, std::shared_ptr<T[]> values_ptr, std::shared_ptr<T[]> b_ptr, std::shared_ptr<T[]> output_dense_ptr);

    // destructor
    virtual ~DMat();

    // print matrix (A contained inside class)
    void PrintMatrix_A();

    // Print vector (b contained inside class)
    void PrintVector_b();

    // Print general
    void PrintVector_x();

    // Call this function if the user has not given an input to generate a random matrix (A) and 
    // vector (b). It takes values for the approx min and max value +- random and a positive definite
    // factor to ensure the random matrix is positive definite.
    // The function can generate symmetric / non-symmetric matrices
    void Populate(T min_val, T max_val, int psd_factor, bool symmetric, int seed=0);

    // solvers
    // Jacobi: Symmetric & Non-Symmetric
    void jacobi(double tolerance=10e-10);

    // G-S: Symmetric (not guaranteed to converge for all Non-Symmetric Matrices)
    void gauss_seidel(double tolerance=10e-10);

    // SOR: Symmetric (not guaranteed to converge for all Non-Symmetric Matrices)
    void successive_over_relaxation(double tolerance=10e-10, double r_factor=1.0, double user_initial_guess=0);

    // Jacobi: Symmetric & Non-Symmetric (explicit)
    void LU_decomp();

    // Jacobi: Symmetric & Non-Symmetric
    void Cholesky();

    // Jacobi: Symmetric & Non-Symmetric
    void Conjugate_Grad(double tolerance=10e-10);

    // methods
    // matrix vector multiplication
    void matA_vec_mult(T *vector_ptr, T *return_vec_ptr);
    T vec_Trans_vec_mult(T *vector_left_ptr, T *vector_right_ptr);

    // The shared pointers to the value, b & output arrays are instantiated here
    // values (matrix A)
    std::shared_ptr<T[]> values = nullptr;
    // b (vector b) (length rows of A)
    std::shared_ptr<T[]> b = nullptr;
    // output_dense (the returned x vector)
    std::shared_ptr<T[]> output_dense = nullptr;

    // rows
    int rows = -1;
    int cols = -1;

    // initialise the min and max values for the random sample matrix 
    // (can be changed by user later)
    T min_val = -10;
    T max_val = 10;
    int psd_factor = 10;

protected:
    bool usr_provided_input_output = false;
    bool symmetric = true;

private:
    // size of values
    int no_values = -1;

};



