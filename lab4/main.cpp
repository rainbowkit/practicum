#include <iostream>
#include <cmath>


double func1(double x) { return std::sin(x); }
double integral1(double x) { return -std::cos(x); }

double (*func)(double)= func1;
double (*integral)(double) = integral1;

double defIntegral(double a, double b) { return integral(b) - integral(a); }

int main() {
    double a, b;
    std::cout << "Enter boundaries: ";
    std::cin >> a >> b;
    while (a > b) {
        std::cout << "Invalid interval, try again: ";
        std::cin >> a >> b;
    }

    double actual = defIntegral(a, b);
    std::cout << "Actual value: " << actual << std::endl;


    return 0;
}
