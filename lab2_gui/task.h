#ifndef TASK_H
#define TASK_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>

double func(double x);

// Класс, позволяющий создать функцию сравнения двух пар (точек) по их удаленности от
// точки x0, при этом без передачи x0 как параметра функции сравнения
// Это нужно чтобы функция сравнения соответствовала формату для std::sort()
// то есть мы можем передать ровно 2 параметра - 2 объекта, которые сравниваются
class CloserTo {
public:
    CloserTo(double x0) {
        this->x0 = x0;
    }
    bool operator()(const std::pair<double, double> &p1, const std::pair<double, double> &p2)
    {
        return std::abs(p1.first - x0) < std::abs(p2.first - x0);
    }

private:
    double x0;
};

/// @brief Генерация таблицы со значениями с фиксированным шагом
/// inverse "переворачивает" таблицу для обратного интерполирования
std::vector<std::pair<double, double>> generateTable(double a, double b, uint16_t nodeCount, bool inverse, double (*func)(double x));

/// @brief Calculate interpolation polynomial in Lagrange form for value x
/// func is used to calculate error
std::pair<double, double> lagrange(double x,
                                   const std::vector<std::pair<double, double>> &table,
                                   uint16_t power,
                                   double (*func)(double x),
                                   bool inverse = false);

/// @brief Calculate interpolation polynomial in Newton form for value x
/// func is used to calculate error
std::pair<double, double> newton(double x,
                                 const std::vector<std::pair<double, double>> &table,
                                 uint16_t power,
                                 double (*func)(double x),
                                 bool inverse = false);

/// @brief Check if vector is sorted by first element of pair
bool isMonot(const std::vector<std::pair<double, double>> &table);

/// @brief find first root of function func(x) = y on interval [a, b] with given error
std::vector<std::pair<double, double>> findFirstRoot(
    const double a,
    const double b,
    uint16_t power,
    const std::vector<std::pair<double, double>> &table,
    double (*func)(double),
    const double y,
    const double error);

#endif // TASK_H
