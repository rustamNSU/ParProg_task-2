#include <iostream>
#include <fstream>
#include<chrono>
#include <omp.h>

#include "BICGStab.hpp"
#include <cmath>

using MatrixType = std::vector<std::vector<double>>;
using VectorType = std::vector<double>;

VectorType matrix_product(const MatrixType &A, const VectorType &x){
    VectorType result(x.size(), 0.0);
    int rowsNumber = A.size();
    int columnsNumber = A[0].size();

    #pragma omp parallel
    {
        #pragma omp for
        for (size_t i = 0; i < rowsNumber; ++i){
            for (size_t j = 0; j < columnsNumber; ++j){
                result[i] += A[i][j] * x[j];
            }
        }
    }
    return result;
}

MatrixType create_test_matrix(int N){
    MatrixType A(N, VectorType(N));
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < N; ++j){
            A[i][j] = 1.0 / ((i-j)*(i-j) - 0.25);
        }
    }
    return A;
}

int main()
{
    using Time = std::chrono::high_resolution_clock;
    using dsec = std::chrono::duration<double>;

    const int g_number_of_threads = 8;
    const int N = 3000;
    omp_set_num_threads(g_number_of_threads);


    MatrixType A = create_test_matrix(N);
    VectorType solution(N, 1.0);
    VectorType rhs = matrix_product(A, solution);

    auto matrixProduct = [&A](const std::vector<double> &x){return matrix_product(A, x);};
    BICGStab<> solver(matrixProduct);



    auto start = Time::now();
    VectorType result = solver.solve(rhs);
    dsec runtime = Time::now() - start;
    auto res = solution - result;
    double error = DotProduct(res, res) / N;


    std::cout << "--- OpenMP in matrix product for BiCGStab ---\n"
              << "    Number of threads    = " << omp_get_max_threads() << '\n'
              << "    Number of processors = " << omp_get_num_procs() << '\n'
              << "    Vector size          = " << N << '\n'
              << "    BiCGStab iterations  = " << solver.get_iterations() << '\n'
              << "    BiCGStab rel. error  = " << solver.get_error() << '\n'
              << "    error (l2-norm)      = " << error << '\n'
              << "    BiCGStab time        = " << runtime.count() << " s.\n"
              << "---------------------------------------------\n";

    std::string filename = "output/openMP_thr_" + std::to_string(omp_get_max_threads()) +
                           "_N_" + std::to_string(N) + ".txt";
    std::ofstream out(filename);
    out << "--- OpenMP in matrix product for BiCGStab ---\n"
        << "    Number of threads    = " << omp_get_max_threads() << '\n'
        << "    Number of processors = " << omp_get_num_procs() << '\n'
        << "    Vector size          = " << N << '\n'
        << "    BiCGStab iterations  = " << solver.get_iterations() << '\n'
        << "    BiCGStab rel. error  = " << solver.get_error() << '\n'
        << "    error (l2-norm)      = " << error << '\n'
        << "    BiCGStab time        = " << runtime.count() << " s.\n"
        << "---------------------------------------------\n";
    out.close();
    return 0;
}
