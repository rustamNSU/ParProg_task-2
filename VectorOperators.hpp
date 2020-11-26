#ifndef TASK_2_VECTOROPERATORS_HPP
#define TASK_2_VECTOROPERATORS_HPP

#include <vector>
#include <cassert>
std::vector<double> operator+(const std::vector<double> &v1, const std::vector<double> &v2)
{
    assert(v1.size()==v2.size() && "Different vectors size");
    auto result = v1;
    for (int i = 0; i < v1.size(); ++i)
    {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

std::vector<double> operator-(const std::vector<double> &v1, const std::vector<double> &v2)
{
    assert(v1.size()==v2.size() && "Different vectors size");
    auto result = v1;
    for (int i = 0; i < v1.size(); ++i)
    {
        result[i] = v1[i] - v2[i];
    }
    return result;
}

std::vector<double> operator*(double a, const std::vector<double> &v)
{
    auto result = v;
    for (int i = 0; i < v.size(); ++i)
    {
        result[i] = a * v[i];
    }
    return result;
}

std::vector<double> operator*(const std::vector<double> &v, double a)
{
    auto result = v;
    for (int i = 0; i < v.size(); ++i)
    {
        result[i] = a * v[i];
    }
    return result;
}

double DotProduct(const std::vector<double> &v1, const std::vector<double> &v2)
{
    assert(v1.size()==v2.size() && "Different vectors size");
    double result = 0.0;
    for (int i = 0; i < v1.size(); ++i)
    {
        result += v1[i] * v2[i];
    }
    return result;
}

#endif //TASK_2_VECTOROPERATORS_HPP
