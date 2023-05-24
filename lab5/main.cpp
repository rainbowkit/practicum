#include <iostream>
#include <boost/math/quadrature/trapezoidal.hpp>
using boost::math::quadrature::trapezoidal;

double func  (const double x) { return sin(x); }
double weight(const double x) { return std::sqrt(x); }
auto f = [](const double x) { return weight(x) * func(x); };

int main() {
    double a, b;
    std::cout << "Enter boundaries: ";
    std::cin >> a >> b;
    while (a > b) {
        std::cout << "Invalid interval, try again: ";
        std::cin >> a >> b;
    }

    std::cout << std::scientific << std::setprecision(12);
    std::cout << "Actual value: " << trapezoidal(f, a, b);


    return 0;
}
