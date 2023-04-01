#include "task.h"

double func(double x) {return sin(x) - x * x / 2;}

/// @brief Генерация таблицы со значениями с фиксированным шагом
std::vector<std::pair<double, double>> generateTable(double a, double b, uint16_t nodeCount) {

    std::vector<std::pair<double, double>> table;
    double step = (b - a) / (nodeCount - 1);
    double x, y;

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Table" << std::endl;

    for (uint8_t i = 0; i < nodeCount - 1; ++i) {
        x = a + i * step;
        y = func(x);
        table.push_back(std::make_pair(x, y));
        std::cout << x << " " << y << std::endl;
    }

    x = b;
    y = func(x);
    table.push_back(std::make_pair(x, y));
    std::cout << x << " " << y << std::endl;

    return table;
}

std::pair<double, double> lagrange(double x, std::vector<std::pair<double, double>> &table, uint16_t power) {
    double sum = 0, summand;
    for (uint16_t k = 0; k <= power; ++k) {
        summand = table.at(k).second;  // f(Xk)
        for (uint16_t i = 0; i <= power; ++i) {
            if (i == k) continue;
            summand *= (x - table.at(i).first) / (table.at(k).first - table.at(i).first);
        }
        sum += summand;
    }
    return std::make_pair(sum, abs(sum - func(x)));
}

std::pair<double, double> newton(double x, std::vector<std::pair<double, double>> &table, uint16_t power) {
    // Таблица разделенных разностей
    std::vector<std::vector<double>> diffTable(table.size() - 1);

    // Разности 1-го порядка
    diffTable.at(0) = std::vector<double>(diffTable.size());
    for (uint16_t i = 0; i < diffTable.at(0).size(); ++i) {
        diffTable.at(0).at(i) = (table.at(i + 1).second - table.at(i).second)/(table.at(i + 1).first - table.at(i).first);
    }

    // Разности порядка > 1
    for (uint16_t i = 1; i < diffTable.size(); ++i) {
        diffTable.at(i) = std::vector<double>(diffTable.size() - i);
        for (uint16_t j = 0; j < diffTable.at(i).size(); ++j) {
            diffTable.at(i).at(j) = (diffTable.at(i - 1).at(j + 1) - diffTable.at(i - 1).at(j)) / (table.at(j + 1 + i).first - table.at(j).first);
        }
    }

    // Вычисление по схеме Горнера
    double sum = diffTable.at(power - 1).at(0);
    for (uint16_t k = power - 1; k > 0; --k) {
        sum *= x - table.at(k).first;
        sum += diffTable.at(k - 1).at(0);
    }
    sum *= x - table.at(0).first;
    sum += table.at(0).second;

    return std::make_pair(sum, abs(sum - func(x)));
}
