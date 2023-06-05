#include <iostream>
#include <cmath>
#include <iomanip>


typedef double(*func_t)(double);
typedef double(*int_func_t)(double, double, double, func_t);

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
double lRect(double a, double b, double steps, func_t f) {
    double step = (b - a) / steps;
    double sum = 0;
    for (uint32_t i = 0; i < steps; ++i) {
        sum += f(a + i * step);
    }
    return step * sum;
}
double rRect(double a, double b, double steps, func_t f) {
    double step = (b - a) / steps;
    double sum = 0;
    for (uint32_t i = 1; i <= steps; ++i) {
        sum += f(a + i * step);
    }
    return step * sum;
}
double mRect(double a, double b, double steps, func_t f) {
    double step = (b - a) / steps;
    double sum = 0;
    for (uint32_t i = 0; i < steps; ++i) {
        sum += f(a + (i + 0.5) * step);
    }
    return step * sum;
}
double trapez(double a, double b, double steps, func_t f) {
    double step = (b - a) / steps;
    double sum = 0;
    for (uint32_t i = 1; i < steps; ++i) {
        sum += f(a + i * step);
    }
    return step * (2 * sum + f(a) + f(b)) / 2;
}
double simpson(double a, double b, double steps, func_t f) {
    double step = (b - a) / steps;
    double sum1 = f(a + step / 2);  // Will be multiplied by 4
    double sum2 = 0;  // Will be multiplied by 2
    for (uint32_t i = 1; i < steps; ++i) {
        sum1 += f(a + step * (0.5 + i));
        sum2 += f(a + step * i);
    }
    return step/6 * (f(a) + f(b) + 4*sum1 + 2*sum2);
}
double int3by8(double a, double b, double steps, func_t f) {
    double step = (b - a) / steps;
    double h = step / 3;
    double sum = 0;
    for (uint32_t i = 0; i < steps; ++i) {
        sum += step * (f(a + i*step) + 3 * (f(a + i*step + h) + f(a + i*step + 2*h)) + f(a + step * (i + 1))) / 8;
    }
    return sum;
}
double runge(uint32_t l, uint32_t d, double int_m, double int_ml) {
    return (std::pow(l,d+1) * int_ml - int_m) / (std::pow(l,d+1) - 1);
}


int main() {
    int_func_t int_func[6] = {lRect, rRect, mRect, trapez, simpson, int3by8};
    std::string func_names[6] = {"Left rect", "Right rect", "Middle rect", "Trapezoidal", "Simpson", "3/8 method"};
    uint8_t precision[6] = {0, 0, 1, 1, 3, 3};

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

    uint32_t l;
    std::cout << "Enter l: ";
    std::cin >> l;

    std::cout << std::scientific << std::setprecision(12);

    double actual = defIntegral(a, b);
    std::cout << "Actual value:\t" << actual << std::endl;
    

    std::cout << "\n\n==========Results===========\n";

    double result1, result2, refined;
    for (uint8_t i = 0; i < 6; ++i) {
        result1 = int_func[i](a, b, steps, func);
        std::cout << func_names[i] << ":\n";
        std::cout << steps << " steps:\t" << result1 << ", Error = " << std::abs(result1 - actual) << std::endl;

        result2 = int_func[i](a, b, steps * l, func);
        std::cout << steps * l << " steps:\t" << result2 << ", Error = " << std::abs(result2 - actual) << std::endl;

        refined = runge(l, precision[i], result1, result2);
        std::cout << "Refined value (Runge):\t" << refined << ", Error = " << std::abs(refined - actual) << "\n\n";
    }

    return 0;
}
