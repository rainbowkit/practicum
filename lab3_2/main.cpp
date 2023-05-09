#include <iostream>
#include <cmath>
#include <vector>

double func(double x) {return     std::exp(3*x);}
double der1(double x) {return 3 * std::exp(3*x);}
double der2(double x) {return 9 * std::exp(3*x);}

// Derivatives for the first and the last points, and any point in between
double derFirst(double y0, double y1, double y2, double step) { return (-3 * y0 + 4 * y1 - y2)/(2 * step); }
double derMiddle(double y0, double y2, double step) { return (y2 - y0)/(2 * step); }
double derLast(double y0, double y1, double y2, double step) { return (3 * y2 - 4 * y1 + y0)/(2 * step); }
double der2Middle(double y0, double y1, double y2, double step) { return (y0 - 2 * y1 + y2)/(step * step); }

std::vector<double> calcDerivatives(const std::vector<double> &x, const std::vector<double> &y) {
    std::vector<double> derivatives;
    double step = x.at(1) - x.at(0);

    derivatives.push_back(derFirst(
            y.at(0),
            y.at(1),
            y.at(2), step));

    for (uint16_t i = 1; i < x.size() - 1; ++i) {
        derivatives.push_back(derMiddle(
                y.at(i - 1),
                y.at(i + 1), step));
    }

    derivatives.push_back(derLast(
            y.at(y.size() - 3),
            y.at(y.size() - 2),
            y.at(y.size() - 1), step));

    return derivatives;
}

std::vector<double> calcDerivatives2(const std::vector<double> &x, const std::vector<double> &y) {
    std::vector<double> derivatives;
    double step = x.at(1) - x.at(0);

    for (uint16_t i = 1; i < x.size() - 1; ++i) {
        derivatives.push_back(der2Middle(
                y.at(i - 1),
                y.at(i),
                y.at(i + 1), step));
    }
    return derivatives;
}

std::pair<std::vector<double>, std::vector<double>> genTable() {
    uint16_t pointCount;
    double a, step;

    std::cout << "Enter number of points: ";
    std::cin >> pointCount;
    while (pointCount < 3) {
        std::cout << "Number of points >= 3 required for this task, try again: ";
        std::cin >> pointCount;
    }

    std::cout << "Enter left bound: ";
    std::cin >> a;

    std::cout << "Enter step: ";
    std::cin >> step;

    std::vector<double> x, y;
    std::cout << "Table of points:\n";
    for (uint16_t i = 0; i < pointCount; ++i) {
        x.push_back(a + i * step);
        y.push_back(func(x.at(i)));
        std::cout << x.at(i) << ": " << y.at(i) << std::endl;
    }
    return std::make_pair(x, y);
}

int main() {
    std::string answer;
    do {
        auto points = genTable();
        auto x = points.first, y = points.second;
        auto derivatives1 = calcDerivatives(x, y);
        auto derivatives2 = calcDerivatives2(x, y);

        std::cout << "x\t\tf(x)\t\tf'(x)\t\tError\t\tf''(x)\t\tError\n" << std::scientific;
        for (uint16_t i = 0; i < x.size(); ++i) {
            std::cout << x.at(i) << '\t'
                      << y.at(i) << '\t'
                      << derivatives1.at(i) << '\t'
                      << std::abs(derivatives1.at(i) - der1(x.at(i))) << '\t';
            if (i > 0 && i < x.size() - 1) {
                std::cout << derivatives2.at(i - 1) << '\t' << std::abs(derivatives2.at(i - 1) - der2(x.at(i)));
            }
            std::cout << std::endl;
        }

        std::cout << "Would you like to enter new parameters? (Yes|no): ";
        std::cin >> answer;
    } while (answer == "yes" || answer == "Yes");


    return 0;
}
