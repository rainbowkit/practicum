#ifndef TASK_H
#define TASK_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>

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
    bool operator() (std::pair<double, double> &p1, std::pair<double, double> &p2) {return abs(p1.first - x0) < abs(p2.first - x0);}
private:
    double x0;
};

/// @brief Генерация таблицы со значениями с фиксированным шагом
std::vector<std::pair<double, double>> generateTable(double a, double b, uint16_t nodeCount);

void calculateValues(std::vector<std::pair<double, double>> &table);
std::pair<double, double> lagrange(double x, std::vector<std::pair<double, double>> &table, uint16_t power);
std::pair<double, double> newton(double x, std::vector<std::pair<double, double>> &table, uint16_t power);

#endif // TASK_H
