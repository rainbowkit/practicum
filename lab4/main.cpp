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

double func7(double x) { return 1.27 * std::pow(x, 5) + 2.04 * x; }
double integral7(double x) { return 1.27/6 * std::pow(x, 6) + 1.02 * std::pow(x, 2); }

double (*func)(double)     =     func7;
double (*integral)(double) = integral7;

double defIntegral(double a, double b) { return integral(b) - integral(a); }


// Quadrature formulas
double lRect(double a, double b, double steps, double(*f)(double)) {
    double step = (b - a) / steps;
    double sum = 0;
    for (uint32_t i = 0; i < steps; ++i) {
        sum += f(a + i * step);
    }
    return step * sum;
}
double rRect(double a, double b, double steps, double(*f)(double)) {
    double step = (b - a) / steps;
    double sum = 0;
    for (uint32_t i = 1; i <= steps; ++i) {
        sum += f(a + i * step);
    }
    return step * sum;
}
double mRect(double a, double b, double steps, double(*f)(double)) {
    double step = (b - a) / steps;
    double sum = 0;
    for (uint32_t i = 0; i < steps; ++i) {
        sum += f(a + (i + 0.5) * step);
    }
    return step * sum;
}
double trapez(double a, double b, double steps, double(*f)(double)) {
    double step = (b - a) / steps;
    double sum = 0;
    for (uint32_t i = 1; i < steps; ++i) {
        sum += f(a + i * step);
    }
    return step * (2 * sum + f(a) + f(b)) / 2;
}
double simpson(double a, double b, double steps, double(*f)(double)) {
    double step = (b - a) / steps;
    double sum1 = f(a + step / 2);  // Will be multiplied by 4
    double sum2 = 0;  // Will be multiplied by 2
    for (uint32_t i = 1; i < steps; ++i) {
        sum1 += f(a + step * (0.5 + i));
        sum2 += f(a + step * i);
    }
    return step/6 * (f(a) + f(b) + 4*sum1 + 2*sum2);
}
double int3by8(double a, double b, double steps, double(*f)(double)) {
    double step = (b - a) / steps;
    double h = step / 3;
    double sum = 0;
    for (uint32_t i = 0; i < steps; ++i) {
        sum += step * (f(a + i*step) + 3 * (f(a + i*step + h) + f(a + i*step + 2*h)) + f(a + step * (i + 1))) / 8;
    }
    return sum;
}


int main() {
    double a, b;
    std::cout << "Enter boundaries: ";
    std::cin >> a >> b;
    while (a > b) {
        std::cout << "Invalid interval, try again: ";
        std::cin >> a >> b;
    }

    uint32_t steps;
    std::cout << "Enter number of steps: ";
    std::cin >> steps;

    double actual = defIntegral(a, b);
    std::cout << "Actual value:\t" << actual << std::endl;

    std::cout << std::scientific;

    double result = lRect(a, b, steps, func);
    std::cout << "Left Rect:\t" << result << ", Error = " << std::abs(result - actual) << std::endl;

    result = rRect(a, b, steps, func);
    std::cout << "Right Rect:\t" << result << ", Error = " << std::abs(result - actual) << std::endl;

    result = mRect(a, b, steps, func);
    std::cout << "Middle Rect:\t" << result << ", Error = " << std::abs(result - actual) << std::endl;

    result = trapez(a, b, steps, func);
    std::cout << "Trapezoidal:\t" << result << ", Error = " << std::abs(result - actual) << std::endl;

    result = simpson(a, b, steps, func);
    std::cout << "Simpson:\t" << result << ", Error = " << std::abs(result - actual) << std::endl;

    result = int3by8(a, b, steps, func);
    std::cout << "3/8 method:\t" << result << ", Error = " << std::abs(result - actual) << std::endl;

    return 0;
}
