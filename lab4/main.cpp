#include <iostream>
#include <cmath>

// Precalculated
double func1(double x) { return std::sin(x); }
double integral1(double x) { return -std::cos(x); }

double func2(double x) { return std::exp(x); }
double (*integral2)(double) = func2;

double func3(double x) { return 3; }
double integral3(double x) { return 3 * x; }

double func4(double x) { return 8 * x; }
double integral4(double x) { return 4 * x * x; }

double func5(double x) { return 9 * x * x; }
double integral5(double x) { return 3 * std::pow(x, 3); }

double func6(double x) { return 4 * std::pow(x, 3); }
double integral6(double x) { return std::pow(x, 4); }

double (*func)(double)     =     func6;
double (*integral)(double) = integral6;

double defIntegral(double a, double b) { return integral(b) - integral(a); }


// Quadrature formulas
double lRect(double a, double b, double(*f)(double)) { return (b - a) * f(a); }
double rRect(double a, double b, double(*f)(double)) { return (b - a) * f(b); }
double mRect(double a, double b, double(*f)(double)) { return (b - a) * f((a + b) / 2) ; }
double trapez(double a, double b, double(*f)(double)) { return (b - a) * (f(a) + f(b)) / 2; }
double simpson(double a, double b, double(*f)(double)) { return (b - a) * (f(a) + 4 * f((a + b) / 2) + f(b)) / 6; }
double int3by8(double a, double b, double(*f)(double)) {
    double h = (b - a) / 3;
    return (b - a) * (f(a) + 3 * (f(a + h) + f(a + 2*h)) + f(b)) / 8;
}


int main() {
    double a, b;
    std::cout << "Enter boundaries: ";
    std::cin >> a >> b;
    while (a > b) {
        std::cout << "Invalid interval, try again: ";
        std::cin >> a >> b;
    }

    double actual = defIntegral(a, b);
    std::cout << "Actual value:\t" << actual << std::endl;

    std::cout << std::scientific;

    double result = lRect(a, b, func);
    std::cout << "Left Rect:\t" << result << ", Error = " << std::abs(result - actual) << std::endl;

    result = rRect(a, b, func);
    std::cout << "Right Rect:\t" << result << ", Error = " << std::abs(result - actual) << std::endl;

    result = mRect(a, b, func);
    std::cout << "Middle Rect:\t" << result << ", Error = " << std::abs(result - actual) << std::endl;

    result = trapez(a, b, func);
    std::cout << "Trapezoidal:\t" << result << ", Error = " << std::abs(result - actual) << std::endl;

    result = simpson(a, b, func);
    std::cout << "Simpson:\t" << result << ", Error = " << std::abs(result - actual) << std::endl;

    result = int3by8(a, b, func);
    std::cout << "3/8 method:\t" << result << ", Error = " << std::abs(result - actual) << std::endl;

    return 0;
}
