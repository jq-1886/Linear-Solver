#include <iostream>
#include <math.h>
#include <ctime>
#include <memory>
#include "DMat.h"
#include "DMat.cpp"
#include <chrono>
#include <fstream>

// This main file tests only the dense methods

// Testing if Ax = b
std::string TestArrays(double *array1, double *array2, int len)
{
    for (int i=0; i<len; i++)
    {
        if (round(array1[i]) != round(array2[i]))
            return "Fail";
    }
    return "Pass";
}

int main()
{
    int loop = 500;
    std::string usr_ans;
    std::cout << "\n \nRunning this will flush out the terminal with a wall of text";
    std::cout << "and may take a minute. Proceed? (y/n) : ";
    std::cin >> usr_ans;
    if (usr_ans == "y")
    {

        int rows = 100;
        int cols = 100;

        double *test_vec = new double [rows] {0};

        int counter = 0; // number of times recorded
        double Jtime = 0; // total time for Jacobi
        double GStime = 0; // total time for Gauss Seidel
        double SORtime = 0; // total time for Successive over relaxation
        double LUtime = 0; // total time for LU decomposition
        double Ctime = 0; // total time for Cholesky
        double CGtime = 0; // total time for Conjugate gradient

        // The average time spent on a method would then be  = the total time recorded / counter

        std::chrono::duration<double> elapsed;
        auto start = std::chrono::high_resolution_clock::now();
        auto finish = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < loop; i++)
        {
            // Instantiating new dense matrix
            auto *dense_sym_matrix = new DMat<double>(rows, cols, false);
            // Randomly populating the matrix with random numbers (using a seed number i
            // to ensure a different matrix every time)
            dense_sym_matrix->Populate(-10, 10, 1000, true, i);


            // The structure to time of the following functions are identical to above 


            // Timing Cholesky
            // Starting timer
            start = std::chrono::high_resolution_clock::now();
            // Apply method
            dense_sym_matrix->Cholesky();
            // End timer
            finish = std::chrono::high_resolution_clock::now();
            elapsed = finish - start;
            dense_sym_matrix->matA_vec_mult(dense_sym_matrix->output_dense.get(), test_vec);
            // If the the input matrix satisfy the conditions by Cholesky, it this matrix will be
            // tested by other functions as well to keep the input of methods consistent.
            // This is required as it is possible some of the random generated matrices are not
            // positive definite, resulting to the methods not function properly and the times recorded
            // will not be reflective of their actual abilities
            if (TestArrays(test_vec, dense_sym_matrix->b.get(), rows) == "Pass")
            {
                std::cout << "\nCholesky : ";

                counter += 1;
                Ctime += elapsed.count();

                // The structure to time of the following functions are identical to above 

                std::cout << "\nLLU decomposition : ";

                // Timing LU decomposition
                start = std::chrono::high_resolution_clock::now();
                dense_sym_matrix->LU_decomp();
                finish = std::chrono::high_resolution_clock::now();
                elapsed = finish - start;
                
                LUtime += elapsed.count();

                std::cout << "\nJacobi : ";

                // Timing Jacobi method
                start = std::chrono::high_resolution_clock::now();
                dense_sym_matrix->jacobi(1e-3);
                finish = std::chrono::high_resolution_clock::now();
                elapsed = finish - start;
                
                Jtime += elapsed.count(); 
                
                std::cout << "\nGauss Seidel : ";

                // Timing Gauss Seidel
                start = std::chrono::high_resolution_clock::now();
                dense_sym_matrix->gauss_seidel(1e-3);
                finish = std::chrono::high_resolution_clock::now();
                elapsed = finish - start;
                
                GStime += elapsed.count();
            
                std::cout << "\nConjugate Gradient : ";

                // Timing Successive over relacation
                start = std::chrono::high_resolution_clock::now();
                dense_sym_matrix->successive_over_relaxation();
                finish = std::chrono::high_resolution_clock::now();
                elapsed = finish - start;

                SORtime += elapsed.count();

                std::cout << "\nConjugate Gradient : ";

                // Timing Conjugate gradient
                start = std::chrono::high_resolution_clock::now();
                dense_sym_matrix->Conjugate_Grad(1e-3);
                finish = std::chrono::high_resolution_clock::now();
                elapsed = finish - start;
                
                CGtime += elapsed.count();
            
            }

            delete dense_sym_matrix;
        }
        
        // Printing out results
        std::cout << std::endl << std::endl;
        std::cout << std::endl << "Averaged over " << counter << " times (out of " << loop << " loops)\n";
        std::cout << "Tolerance = " << 1e-3 << std::endl;
        std::cout << "Input Matrix A size is " << rows << " x " << rows << std::endl;
        std::cout << "\nAverage time for Jacobi = " << Jtime/counter*1000 << "milliseconds";
        std::cout << "\nAverage time for Gauss Seidel = " << GStime/counter*1000 << "milliseconds";
        std::cout << "\nAverage time for Successive over relaxation = " << SORtime/counter*1000 << "milliseconds";
        std::cout << "\nAverage time for Conjugate Gradient = " << CGtime/counter*1000 << "milliseconds";
        std::cout << "\nAverage time for LU = " << LUtime/counter*1000 << "milliseconds";
        std::cout << "\nAverage time for Cholesky = " << Ctime/counter*1000 << "milliseconds"<< std::endl;

        // Putting results into text file: timing_results.txt
        std::ofstream textfile;
        textfile.open("timing_results.txt");
        textfile << std::endl << "Averaged over " << counter << " times (out of " << loop << " loops)\n";
        textfile << "Tolerance = " << 1e-3 << std::endl;
        textfile << "Input Matrix A size is " << rows << " x " << rows << std::endl;
        textfile << "\nAverage time for Jacobi = " << Jtime/counter*1000 << "milliseconds";
        textfile << "\nAverage time for Gauss Seidel = " << GStime/counter*1000 << "milliseconds";
        textfile << "\nAverage time for Successive over relaxation = " << SORtime/counter*1000 << "milliseconds";
        textfile << "\nAverage time for Conjugate Gradient = " << CGtime/counter*1000 << "milliseconds";
        textfile << "\nAverage time for LU = " << LUtime/counter*1000 << "milliseconds";
        textfile << "\nAverage time for Cholesky = " << Ctime/counter*1000 << "milliseconds"<< std::endl;

    }
    else
    {
        return 0;
    }
    return 0;
}