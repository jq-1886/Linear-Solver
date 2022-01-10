#include <iostream>
#include <math.h>
#include <ctime>
#include <memory>
#include "CSRMat.h"
#include "CSRMat.cpp"
#include "DMat.h"
#include "DMat.cpp"
#include <vector>

// Function to test if two arrays are the same (of templated type)
template <class T>
std::string TestArrays(T *array1, T *array2, int len)
{
    for (int i=0; i<len; i++)
    {
        if (abs(array1[i] - array2[i]) > 10e-6)
            return "Fail";
    }
    return "Pass";
}

// Function to test if two vectors are the same (of templated type)
template <class T>
std::string TestVectors(std::vector<T>* vec1, std::vector<T>* vec2, int len)
{
    for (int i=0; i< len; i++)
    {
        if (abs(vec1->at(i) - vec2->at(i)) > 10e-4)
            return "Fail";
    }
    return "Pass";
}

int main()
{
    // First look at the Dense matrix solver Class
    // we will automatically generate a a random matrix of 15 x 15

    int rows = 15;
    int cols = 15;

    // Instantiate the Dense matrix class and specifiy no user input (last parameter false)
    // Use a shared pointer when instantiating the class to avoid having to delete it at the end
    auto dense_sym_matrix = std::make_shared<DMat<double>>(rows, cols, false);
    // Call the populate function to generate the matrix A and vector b
    // Generate a symmetric matrix (last parameter = true)
    dense_sym_matrix->Populate(-10, 10, 100, true);

    // ---------------------------------------------------------------------------
    // If the user wishes to input a matrix then call the class as follows:
    // auto dense_user_matrix = std::make_shared<DMat<double>>(rows, cols, user_values_ptr, user_RHS_ptr, user_output_ptr);
    // where the pointers are shared pointers.
    // ---------------------------------------------------------------------------

    // Option to print the generated matrices
    // dense_sym_matrix->PrintMatrix_A();
    // dense_sym_matrix->PrintVector_b();


    // Method to Call Solvers for Real Positive Definite Matrices
    // we are creating a test vector which will be the result of A . x which is then compared against
    // the input b vector to see if the solution is accurate
    double *test_vec = new double [rows] {0};

    // Solvers for Symmetric Dense Matrices:
    std::cerr << std::endl << "\t Dense solvers on Randomly Generated Symmetric Dense Matrix size: " << rows << ", " << cols << std::endl << std::endl;

    std::cerr << std::endl << "Jacobi Solver: ";
    dense_sym_matrix->jacobi();
    dense_sym_matrix->matA_vec_mult(dense_sym_matrix->output_dense.get(), test_vec);
    std::cerr << "\t" << TestArrays(test_vec, dense_sym_matrix->b.get(), rows) << std::endl;
    dense_sym_matrix->PrintVector_x();

    std::cerr << std::endl << "GS Solver: ";
    dense_sym_matrix->gauss_seidel();
    dense_sym_matrix->matA_vec_mult(dense_sym_matrix->output_dense.get(), test_vec);
    std::cerr << "\t" << TestArrays(test_vec, dense_sym_matrix->b.get(), rows) << std::endl;
    dense_sym_matrix->PrintVector_x();

    std::cerr << std::endl << "SOR Solver: ";
    dense_sym_matrix->successive_over_relaxation(10e-10);
    dense_sym_matrix->matA_vec_mult(dense_sym_matrix->output_dense.get(), test_vec);
    std::cerr << "\t" << TestArrays(test_vec, dense_sym_matrix->b.get(), rows) << std::endl;
    dense_sym_matrix->PrintVector_x();
    
    std::cerr << std::endl << "LU Solver: ";
    dense_sym_matrix->LU_decomp();
    dense_sym_matrix->matA_vec_mult(dense_sym_matrix->output_dense.get(), test_vec);
    std::cerr << "\t" << TestArrays(test_vec, dense_sym_matrix->b.get(), rows) << std::endl;
    dense_sym_matrix->PrintVector_x();

    std::cerr << std::endl << "Conjugate Gradient: ";
    dense_sym_matrix->Conjugate_Grad();
    dense_sym_matrix->matA_vec_mult(dense_sym_matrix->output_dense.get(), test_vec);
    std::cerr << "\t" << TestArrays(test_vec, dense_sym_matrix->b.get(), rows) << std::endl;
    dense_sym_matrix->PrintVector_x();

    std::cerr << std::endl << "Cholesky: ";
    dense_sym_matrix->Cholesky();
    dense_sym_matrix->matA_vec_mult(dense_sym_matrix->output_dense.get(), test_vec);
    std::cerr << "\t" << TestArrays(test_vec, dense_sym_matrix->b.get(), rows) << std::endl;
    dense_sym_matrix->PrintVector_x();
    
    std::cerr << std::endl;


    // Solvers for Non-Symmetric Dense Matrices
    // Instantiate the Dense matrix class and specifiy no user input (last parameter false)
    // Use a shared pointer when instantiating the class to avoid having to delete it at the end
    auto dense_non_sym_matrix = std::make_shared<DMat<double>>(rows, cols, false);
    // Call the populate function to generate the matrix A and vector b
    // Generate a symmetric matrix (last parameter = false)
    dense_non_sym_matrix->Populate(-10, 10, 100, false);

    std::cerr << std::endl << "\t Dense solvers on Randomly Generated Non-Symmetric Dense Matrix size: " << rows << ", " << cols << std::endl << std::endl;

    std::cerr << std::endl << "Jacobi Solver: ";
    dense_non_sym_matrix->jacobi(10e-10);
    dense_non_sym_matrix->matA_vec_mult(dense_non_sym_matrix->output_dense.get(), test_vec);
    std::cerr << "\t" << TestArrays(test_vec, dense_non_sym_matrix->b.get(), rows) << std::endl;
    dense_non_sym_matrix->PrintVector_x();

    std::cerr << std::endl << "Conjugate Gradient: ";
    dense_non_sym_matrix->Conjugate_Grad();
    dense_non_sym_matrix->matA_vec_mult(dense_non_sym_matrix->output_dense.get(), test_vec);
    std::cerr << "\t" << TestArrays(test_vec, dense_non_sym_matrix->b.get(), rows) << std::endl;
    dense_non_sym_matrix->PrintVector_x();

    std::cerr << std::endl << "LU Solver: ";
    dense_non_sym_matrix->LU_decomp();
    dense_non_sym_matrix->matA_vec_mult(dense_non_sym_matrix->output_dense.get(), test_vec);
    std::cerr << "\t" << TestArrays(test_vec, dense_non_sym_matrix->b.get(), rows) << std::endl;
    dense_non_sym_matrix->PrintVector_x();

    delete[] test_vec;

    // ---------------------------------------------------------------------
    // Solvers for Symmetric Sparse Matrices
    // ---------------------------------------------------------------------

    // User defined vectors for the system stored in CSR
    int rows1 = 11;
    int cols1 = 11;
    int nnzs1 = 41;
    std::vector<long double> user_values1 {20, 1, 4, 20, 1, 2, 3, 20, 1, 8, 20, 6, 1, 1, 20, 1, 6, 1, 6, 20, 5, 4, 8, 20, 1,1, 5, 1, 20, 1, 2, 2, 1, 20, 3, 6, 2, 20, 5, 5, 20};
    std::vector<int> user_col_idx1 {0, 4, 6, 1, 4, 8, 9, 2, 5, 6, 3, 5, 0, 1, 4, 7, 9, 2, 3, 5, 7, 0, 2, 6, 7, 4, 5, 6, 7, 8, 9, 1, 7, 8, 1, 4, 7, 9, 10, 9, 10};
    std::vector<int> user_row_pos1 {0, 3, 7, 10, 12, 17, 21, 25, 31, 34, 39, 41};
    std::vector<long double> user_RHS1 {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    std::vector<long double> user_output1 (rows, 0);

    // The user pointers to the vectors
    std::vector<long double>* user_values_ptr1 = &user_values1;
    std::vector<int>* user_col_ptr1 = &user_col_idx1;
    std::vector<int>* user_row_ptr1 = &user_row_pos1;
    std::vector<long double>* user_RHS_ptr1 = &user_RHS1;
    std::vector<long double>* user_output_ptr1 = &user_output1;

    // instantiate the class using the constructor (overloaded)
    auto *sparse_user_matrix = new CSRMat<long double>(rows1, cols1, nnzs1, user_values_ptr1, user_row_ptr1, user_col_ptr1, user_RHS_ptr1, user_output_ptr1);

    ///// Testing framework //////
    // initialise a testing vector to check the product of A . x
    std::vector<long double>* test_vec_ptr1 = new std::vector<long double>(rows1, 0);

    std::cerr << std::endl << "\t Sparse Solvers on User Defined Symmetric Sparse Matrix size: " << rows1 << ", " << cols1 << std::endl << std::endl;

    std::cout << "Jacobi Solver: ";
    sparse_user_matrix->jacobi(10e-20);
    sparse_user_matrix->matA_vec_mult(sparse_user_matrix->output_sparse, test_vec_ptr1);
    std::cerr << "\t" << TestVectors(test_vec_ptr1, sparse_user_matrix->b, rows1) << std::endl;
    sparse_user_matrix->PrintVector_x();

    std::cout << "Gauss Seidel: ";
    sparse_user_matrix->gauss_seidel(10e-10);
    sparse_user_matrix->matA_vec_mult(sparse_user_matrix->output_sparse, test_vec_ptr1);
    std::cerr << "\t" << TestVectors(test_vec_ptr1, sparse_user_matrix->b, rows1) << std::endl;
    sparse_user_matrix->PrintVector_x();

    std::cout << "Successive Over Relaxation: ";
    sparse_user_matrix->successive_over_relaxation(10e-10);
    sparse_user_matrix->matA_vec_mult(sparse_user_matrix->output_sparse, test_vec_ptr1);
    std::cerr << "\t" << TestVectors(test_vec_ptr1, sparse_user_matrix->b, rows1) << std::endl;
    sparse_user_matrix->PrintVector_x();

    std::cout << "Cholesky: ";
    sparse_user_matrix->Cholesky();
    sparse_user_matrix->matA_vec_mult(sparse_user_matrix->output_sparse, test_vec_ptr1);
    std::cerr << "\t" << TestVectors(test_vec_ptr1, sparse_user_matrix->b, rows1) << std::endl;
    sparse_user_matrix->PrintVector_x();

    std::cout << "Conjugate Gradient: ";
    sparse_user_matrix->Conjugate_Grad(10e-10);
    sparse_user_matrix->matA_vec_mult(sparse_user_matrix->output_sparse, test_vec_ptr1);
    std::cerr << "\t" << TestVectors(test_vec_ptr1, sparse_user_matrix->b, rows1) << std::endl;
    sparse_user_matrix->PrintVector_x();

    delete test_vec_ptr1;


    // ---------------------------------------------------------------------
    // Solvers for Non-Symmetric Sparse Matrices
    // ---------------------------------------------------------------------

    // User defined vectors for a Non-Symmetric Sparse CSR format
    int rows2 = 11;
    int cols2 = 11;
    int nnzs2 = 39;
    std::vector<long double> user_values2 {20, 4, 20, 1, 2, 3, 20, 1, 8, 20, 6, 1, 1, 20, 1, 6, 1, 6, 20, 5, 4, 8, 20, 1, 1, 1, 20, 1, 2, 1, 20, 3, 6, 2, 20, 5, 8, 5, 20};
    std::vector<int> user_col_idx2 {0, 5, 1, 4, 8, 9, 2, 5, 6, 3, 5, 0, 1, 4, 7, 9, 2, 3, 5, 7, 0, 2, 6, 7, 4, 6, 7, 8, 9, 7, 8, 1, 4, 7, 9, 10, 3, 9, 10};
    std::vector<int> user_row_pos2 {0, 2, 6, 9, 11, 16, 20, 24, 29, 31, 36, 39};
    std::vector<long double> user_RHS2 {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    std::vector<long double> user_output2 (rows, 0);

    // The user pointers to the vectors
    std::vector<long double>* user_values_ptr2 = &user_values2;
    std::vector<int>* user_col_ptr2 = &user_col_idx2;
    std::vector<int>* user_row_ptr2 = &user_row_pos2;
    std::vector<long double>* user_RHS_ptr2= &user_RHS2;
    std::vector<long double>* user_output_ptr2 = &user_output2;


    // instantiate the class using the constructor (overloaded)
    auto *sparse_user_matrix2 = new CSRMat<long double>(rows2, cols2, nnzs2, user_values_ptr2, user_row_ptr2, user_col_ptr2, user_RHS_ptr2, user_output_ptr2);

    ///// Testing framework //////
    // initialise a testing vector to check the product of A . x
    std::vector<long double>* test_vec_ptr2 = new std::vector<long double>(rows2, 0);

    std::cerr << std::endl << "\t Sparse Solvers on User Defined Non-Symmetric Sparse Matrix size: " << rows2 << ", " << cols2 << std::endl << std::endl;

    // This works for this test cas but in general Jacobi is not guaranteed to converge for Non-Symmetric matrices
    std::cout << "Jacobi Solver: ";
    sparse_user_matrix2->jacobi(10e-20);
    sparse_user_matrix2->matA_vec_mult(sparse_user_matrix2->output_sparse, test_vec_ptr2);
    std::cerr << "\t" << TestVectors(test_vec_ptr2, sparse_user_matrix2->b, rows2) << std::endl;
    sparse_user_matrix2->PrintVector_x();




    delete test_vec_ptr2;

    // ---------------------------------------------------------------------
    // Solvers for Non-Symmetric Banded Sparse Matrices
    // ---------------------------------------------------------------------

    // User defined vectors for a Banded Non-Symmetric Sparse CSR format
    int rows3 = 5;
    int cols3 = 5;
    int nnzs3 = 13;
    std::vector<long double> user_values3 {20, 4, 4, 20, 3, 8, 20, 7, 7, 20, 2, 4, 20};
    std::vector<int> user_col_idx3 {0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4};
    std::vector<int> user_row_pos3 {0, 2, 5, 8, 11, 13};
    std::vector<long double> user_RHS3 {1, 2, 3, 4, 5};
    std::vector<long double> user_output3 (rows, 0);

    // The user pointers to the vectors
    std::vector<long double>* user_values_ptr3 = &user_values3;
    std::vector<int>* user_col_ptr3 = &user_col_idx3;
    std::vector<int>* user_row_ptr3 = &user_row_pos3;
    std::vector<long double>* user_RHS_ptr3 = &user_RHS3;
    std::vector<long double>* user_output_ptr3 = &user_output3;

    // instantiate the class using the constructor (overloaded)
    auto *sparse_user_matrix3 = new CSRMat<long double>(rows3, cols3, nnzs3, user_values_ptr3, user_row_ptr3, user_col_ptr3, user_RHS_ptr3, user_output_ptr3);

    ///// Testing framework //////
    // initialise a testing vector to check the product of A . x
    std::vector<long double>* test_vec_ptr3 = new std::vector<long double>(rows3, 0);

    std::cerr << std::endl << "\t Sparse Solvers on User Defined Banded Non-Symmetric Sparse Matrix size: " << rows3 << ", " << cols3 << std::endl << std::endl;

    std::cout << "LU Banded: ";
    sparse_user_matrix3->LU_Banded();
    sparse_user_matrix3->matA_vec_mult(sparse_user_matrix3->output_sparse, test_vec_ptr3);
    std::cerr << "\t" << TestVectors(test_vec_ptr3, sparse_user_matrix3->b, rows3) << std::endl;
    sparse_user_matrix3->PrintVector_x();

    delete test_vec_ptr3;




    // class shared pointers will fall out of scope


    return 0;
}