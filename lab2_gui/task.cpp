#include "task.h"

double func(double x)
{
    return sin(x) - x * x / 2;
}

/// @brief Генерация таблицы со значениями с фиксированным шагом
std::vector<std::pair<double, double>> generateTable(double a, double b, uint16_t nodeCount, bool inverse, double (*func)(double x)) {

    std::vector<std::pair<double, double>> table;
    double step = (b - a) / (nodeCount - 1);
    double x, y;

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Table" << std::endl;

    for (uint8_t i = 0; i < nodeCount - 1; ++i) {
        x = a + i * step;
        y = func(x);
        table.push_back((inverse ? std::make_pair(y, x) : std::make_pair(x, y)));
        std::cout << x << " " << y << std::endl;
    }

    x = b;
    y = func(x);
    table.push_back((inverse ? std::make_pair(y, x) : std::make_pair(x, y)));
    std::cout << x << " " << y << std::endl;

    return table;
}

std::pair<double, double> lagrange(double x,
                                   const std::vector<std::pair<double, double>> &table,
                                   uint16_t power,
                                   double (*func)(double x),
                                   bool inverse)
{
    double sum = 0, summand;
    for (uint16_t k = 0; k <= power; ++k) {
        summand = table.at(k).second;  // f(Xk)
        for (uint16_t i = 0; i <= power; ++i) {
            if (i == k) continue;
            summand *= (x - table.at(i).first) / (table.at(k).first - table.at(i).first);
        }
        sum += summand;
    }
    return std::make_pair(sum, (inverse ? std::abs(func(sum) - x) : std::abs(sum - func(x))));
}

std::pair<double, double> newton(double x,
                                 const std::vector<std::pair<double, double>> &table,
                                 uint16_t power,
                                 double (*func)(double x),
                                 bool inverse)
{
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

    return std::make_pair(sum, (inverse ? std::abs(func(sum) - x) : std::abs(sum - func(x))));
}

/// @brief Отделение корней способом табулирования
/// @param a Левая граница отрезка
/// @param b Правая граница отрезка
/// @param numOfSteps Количество шагов, отрезков, на которых проверяется перемена знака
/// @param func Указатель на саму функцию
/// @return Список отрезков перемены знака
std::vector<std::pair<double, double>> separateRoots(double a, double b, uint16_t numOfSteps, std::function<double(double)>& func) {
    std::vector<std::pair<double, double>> roots;
    double step = (b - a) / numOfSteps;
    double x1 = a;
    double x2 = a + step;
    double y1 = func(x1), y2 = func(x2);

    std::cout << "Searching for roots in [" << a << ", " << b << "] ..." << std::endl;
    for (uint16_t i = 0; i < numOfSteps; ++i) {

        if (y1 * y2 < 0) {
            roots.push_back(std::make_pair(x1, x2));
            std::cout << "Found root between " << x1 << " and " << x2 << std::endl;
        }
        if (y2 != 0) {  // Если y2 == 0, осталось проверить четность этого
            x1 = x2;    // корня, подвинув только правую границу отрезка
            y1 = y2;
        }
        x2 = a + (i + 2) * step;  // x2 += step;
        y2 = func(x2);
    }
    std::cout << "Number of roots: " << roots.size() << std::endl;
    return roots;
}


/// @brief Метод секущих
/// @param roots Список отрезков перемены знака
/// @param error Точность
/// @param func  Указатель на саму функцию
/// @return Список пар вида (приближённое значение, погрешность)
std::vector<std::pair<double, double>> secantMethod(
    const std::vector<std::pair<double, double>> &roots, double error, std::function<double(double)>& func)
{
    std::vector<std::pair<double, double>> refined_roots;
    double x1, x2, derivative;
    uint16_t steps;
    std::cout << "\n\nSecant:\n=====\n";

    for (uint16_t i = 0; i < roots.size(); ++i) {
        x1 = (roots[i].first + roots[i].second) / 2;
        derivative = (func(roots[i].second) - func(roots[i].first)) / (roots[i].second - roots[i].first);
        x2 = x1 - func(x1) / derivative;
        steps = 1;
        while (std::abs(x2 - x1) > 2 * error) {
            ++steps;

            derivative = (func(x2) - func(x1)) / (x2 - x1);
            x1 = x2;
            x2 = x1 - func(x1) / derivative;
        }
        refined_roots.push_back(std::make_pair(x2, std::abs(x2 - x1) / 2));
        std::cout << "Root " << i << ": " << steps << " steps,\t";
        std::cout << "|Xm - Xm-1| = " << std::abs(x2 - x1);
        std::cout << ", |f(x) - 0| = " << std::abs(func(x2)) << std::endl;
    }

    std::cout << "=====\n";
    return refined_roots;
}

bool isMonot(const std::vector<std::pair<double, double>> &table)
{
    bool acsending = table.at(0).first < table.at(1).first;
    for (uint16_t i = 1; i < table.size() - 1; ++i) {
        if (acsending != (table.at(i).first < table.at(i + 1).first)
            || table.at(i).first == table.at(i + 1).first) {
            return false;
        }
    }
    return true;
}

std::vector<std::pair<double, double>> findFirstRoot(
    const double a,
    const double b,
    uint16_t power,
    const std::vector<std::pair<double, double>> &table,
    double (*func)(double),
    const double y,
    const double error)
{
    std::function<double(double)> shiftFunc = [&func, &table, power, y](const double x) {
        return lagrange(x, table, power, func).first - y;
    };
    auto roots = separateRoots(a, b, 100, shiftFunc);
    if (roots.size() == 0) { return std::vector<std::pair<double, double>>(); }
    return secantMethod(roots, error, shiftFunc);
}
