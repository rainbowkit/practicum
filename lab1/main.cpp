/*
    Никита Ардашев, 21.Б06-мм
    Задание 1
    Численные методы решения нелинейных уравнений
*/

#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

// Функция, для которой осуществляется поиск приближенных значений корней
double func(double x) { return x - 10 * sin(x); }

/// @brief Производная функции func в точке x
double get_derivative(double (*func)(double), double x) {
    // double dx = 0.0000000001;
    // return (func(x + dx) - func(x)) / dx;
    return 1 - 10 * cos(x);
}

/// @brief Отделение корней способом табулирования
/// @param a Левая граница отрезка
/// @param b Правая граница отрезка 
/// @param numOfSteps Количество шагов, отрезков, на которых проверяется перемена знака
/// @param func Указатель на саму функцию
/// @return Список отрезков перемены знака
vector<pair<double, double>> separate_roots(double a, double b, uint16_t numOfSteps, double (*func)(double x)) {
    vector<pair<double, double>> roots;
    double step = (b - a) / numOfSteps;
    double x1 = a;
    double x2 = a + step; 
    double y1 = func(x1), y2 = func(x2);

    cout << "Searching for roots in [" << a << ", " << b << "] ..." << endl;
    for (uint16_t i = 0; i < numOfSteps; ++i) {

        if (y1 * y2 < 0) {
            roots.push_back(make_pair(x1, x2));
            cout << "Found root between " << x1 << " and " << x2 << endl;
        }
        if (y2 != 0) {  // Если y2 == 0, осталось проверить четность этого
            x1 = x2;    // корня, подвинув только правую границу отрезка
            y1 = y2;
        } 
        x2 = a + (i + 2) * step;  // x2 += step;
        y2 = func(x2);
    }
    cout << "Number of roots: " << roots.size() << endl;
    return roots;
}

/// @brief Метод бисекции
/// @param roots Список отрезков перемены знака
/// @param error Точность
/// @param func  Указатель на саму функцию
/// @return Список пар вида (приближенное значение, погрешность)
vector<pair<double, double>> bisection(vector<pair<double, double>> &roots, double error, double (*func)(double)) {
    double a, b, c, y_c;
    uint16_t steps;
    vector<pair<double, double>> refined_roots;
    cout << "\n\nBisection:\n=====\n";

    for (uint16_t i = 0; i < roots.size(); ++i) {
        steps = 0;
        a = roots[i].first;
        b = roots[i].second;
        
        while (b - a > 2 * error) {
            ++steps;

            c = (a + b) / 2;
            y_c = func(c);

            if (func(a) * y_c < 0) { b = c; }
            else if (y_c == 0) {  // Exact root
                a = y_c;
                b = y_c;
                break;
            }
            else { a = c; }
        }

        refined_roots.push_back(make_pair(
            (a + b) / 2,
            (b - a) / 2
        ));
        cout << "Root " << i << ": " << steps << " steps,\t";
        cout << "|Xm - Xm-1| = " << b - a;
        cout << ", |f(x) - 0| = " << abs(func(refined_roots[i].first)) << endl;
    }
    cout << "=====\n";
    return refined_roots;
}

/// @brief Обычный и модифицированный метод Ньютона
/// @param roots Список отрезков перемены знака
/// @param error Точность
/// @param func  Указатель на саму функцию
/// @param modified true - использовать модифицированный метод, false - обычный
/// @return Список пар вида (приближённое значение, погрешность)
vector<pair<double, double>> newton_roots(vector<pair<double, double>> &roots, double error, double (*func)(double), bool modified) {
    vector<pair<double, double>> refined_roots;
    double x1, x2, derivative;
    uint16_t steps;
    cout << "\n\nNewton";
    if (modified) cout << " modified";
    cout << ":\n=====\n";

    for (uint16_t i = 0; i < roots.size(); ++i) {
        x1 = (roots[i].first + roots[i].second) / 2;
        derivative = get_derivative(func, x1);
        x2 = x1 - func(x1) / derivative;
        steps = 1;
        while (abs(x2 - x1) > 2 * error) {
            ++steps;
            
            x1 = x2;
            x2 = x1 - func(x1) / (modified ? derivative : get_derivative(func, x1));
        }
        refined_roots.push_back(make_pair(x2, abs(x2 - x1) / 2));
        cout << "Root " << i << ": " << steps << " steps,\t";
        cout << "|Xm - Xm-1| = " << abs(x2 - x1);
        cout << ", |f(x) - 0| = " << abs(func(x2)) << endl;
    }

    cout << "=====\n";
    return refined_roots;
}

/// @brief Метод секущих
/// @param roots Список отрезков перемены знака
/// @param error Точность
/// @param func  Указатель на саму функцию
/// @return Список пар вида (приближённое значение, погрешность)
vector<pair<double, double>> secant_method(vector<pair<double, double>> &roots, double error, double (*func)(double)) {
    vector<pair<double, double>> refined_roots;
    double x1, x2, derivative;
    uint16_t steps;
    cout << "\n\nSecant:\n=====\n";

    for (uint16_t i = 0; i < roots.size(); ++i) {
        x1 = (roots[i].first + roots[i].second) / 2;
        derivative = (func(roots[i].second) - func(roots[i].first)) / (roots[i].second - roots[i].first);
        x2 = x1 - func(x1) / derivative;
        steps = 1;
        while (abs(x2 - x1) > 2 * error) {
            ++steps;
            
            derivative = (func(x2) - func(x1)) / (x2 - x1);
            x1 = x2;
            x2 = x1 - func(x1) / derivative;
        }
        refined_roots.push_back(make_pair(x2, abs(x2 - x1) / 2));
        cout << "Root " << i << ": " << steps << " steps,\t";
        cout << "|Xm - Xm-1| = " << abs(x2 - x1);
        cout << ", |f(x) - 0| = " << abs(func(x2)) << endl;
    }

    cout << "=====\n";
    return refined_roots;
}

// Выводит в консоль корни и соответствующие погрешности
void print_roots(const vector<pair<double, double>> &roots) {
    cout << "Roots: " << endl;
    for (auto root : roots) {
        if (root.first >= 0) cout << ' ';
        cout << root.first << " +- " << root.second << endl;
    }
}


int main() {
    cout << "=== Task 1 ===\n";
    cout << "Function: x - 10 * sin(x)\n";

    // Границы отрезка, точность и количество шагов
    double a = -10, b = 10, error = pow(10.0, -10);
    uint16_t numOfSteps = 1000;

    // Root separation
    auto roots = separate_roots(a, b, numOfSteps, &func);

    // Root refinement
    cout << fixed << setprecision(17);

    auto refined_roots = bisection(roots, error, &func);
    print_roots(refined_roots);

    refined_roots = newton_roots(roots, error, &func, false);
    print_roots(refined_roots);

    refined_roots = newton_roots(roots, error, &func, true);
    print_roots(refined_roots);

    refined_roots = secant_method(roots, error, &func);
    print_roots(refined_roots);

    return 0;
}