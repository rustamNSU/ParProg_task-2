#ifndef TASK_2_BICGSTAB_HPP
#define TASK_2_BICGSTAB_HPP

#include <functional>
#include <cassert>
#include <vector>
#include <cmath>
#include <limits>
#include <chrono>

#include "VectorOperators.hpp"

template <class VectorType = std::vector<double>>
class BICGStab{
private:
    using VectorFunction = std::function<VectorType(const VectorType&)>;
    VectorFunction matrixProduct_;   // defined A*x (class BICGStab does not know implementation)

    int iterations_ = 0;
    int maxIterations_;
    double error_;
    double tolerance_ = 10*std::numeric_limits<double>::epsilon();

public:
    BICGStab() = delete;
    BICGStab(VectorFunction matrixProduct): matrixProduct_(matrixProduct){

    }
    int get_iterations() const{
        return iterations_;
    }
    double get_error() const{
        return error_;
    }
    VectorType solve(const VectorType &rhs){
        VectorType initial(rhs.size(), 0.0);
        return solve(rhs, initial, tolerance_);
    }
    VectorType solve(const VectorType &rhs, VectorType x, double tolerance){
        int n = rhs.size();
        maxIterations_ = n;

        VectorType r  = rhs - matrixProduct_(x);
        VectorType r0 = r;
        double rho    = 1.0;
        double alpha  = 1.0;
        double omega  = 1.0;

        VectorType v(n, 0.0);
        VectorType p(n, 0.0);
        VectorType s(n, 0.0);
        VectorType t(n, 0.0);

        /* using relative tolerance */
        double epsilon = tolerance * tolerance * DotProduct(rhs, rhs);
        int i = 0;

        while(DotProduct(r, r) > epsilon && i < maxIterations_){
            double rho_old = rho;
            rho = DotProduct(r0, r);
            double beta = (rho / rho_old) * (alpha / omega);
            p = r + beta * (p - omega * v);

            /* A*p */
            v = matrixProduct_(p);
            alpha = rho / DotProduct(r0, v);
            s = r - alpha * v;
            t = matrixProduct_(s);
            double t_sqnorm = DotProduct(t, t);
            if ( t_sqnorm > 0.0 ){
                omega = DotProduct(t, s) / t_sqnorm;
            }
            else {
                omega = 0.0;
            }
            x = x + alpha * p + omega * s;
            r = s - omega * t;
            ++i;
        }
        error_ = std::sqrt(DotProduct(r, r) / DotProduct(rhs, rhs));  // relative error
        iterations_ = i;
        return x;
    }

    ~BICGStab() = default;
};




#endif //TASK_2_BICGSTAB_HPP
