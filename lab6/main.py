import sys
import math
import os
import sympy
from prettytable import PrettyTable

x = sympy.Symbol('x')
y = sympy.Symbol('y')

ITER_MAX = 100
Y_DERIVATIVE = "-y(x)+x"
Y_SOLUTION = "2e^(-x)+x-1"
Y = 2 * sympy.exp(-x) + x - 1
F = -y + x
Y_0 = 1
X_0 = 0
TAYLOR_DEGREE = 6
ADAMS_COEFS = [1, 1 / 2, 5 / 12, 3 / 8, 251 / 720]
ROUNDING = 8


def get_Taylor_multiplier(x_0, n):
    y_deriv = Y
    y_deriv = sympy.diff(y_deriv, x, n)
    return (y_deriv.subs(x, X_0) / math.factorial(n)) * math.pow(x_0 - X_0, n)


def get_Taylor_series(x_0):
    result = 0
    for i in range(TAYLOR_DEGREE):
        result += get_Taylor_multiplier(x_0, i)
    return result


def Taylor_series_method(accurate_values, h, N):
    table = PrettyTable()
    table.field_names = ["k", "x_k", "Приближенное значение y(x_k)", "Абсолютная погрешность"]
    Taylor_values = []
    for k in range(-2, N + 1):
        x_k = accurate_values[k + 2][0]
        y_k = get_Taylor_series(x_k)
        table.add_row([k, x_k, y_k, abs(y_k - accurate_values[k + 2][1])])
        Taylor_values.append([x_k, y_k])
    print("\n\nМетод разложения в ряд Тейлора:")
    print(table)
    return Taylor_values


def get_div_diff(etas, i, j):
    if j == 0:
        return etas[i]
    if j == 1:
        return etas[i - 1] - etas[i]
    return get_div_diff(etas, i - 1, j - 1) - get_div_diff(etas, i, j - 1)


def get_Adams_value(known_values, h):
    result = known_values[-1][1]
    etas = []
    for n in range(5):
        etas.append(h * F.subs({x: known_values[-(n + 1)][0], y: known_values[-(n + 1)][1]}))
    for i in range(5):
        result += ADAMS_COEFS[i] * get_div_diff(etas, i, i)
    return result


def Adams_method(accurate_values, Taylor_values, h, N):
    table = PrettyTable()
    table.field_names = ["k", "x_k", "Приближенное значение y(x_k)", "Абсолютная погрешность"]
    Adams_values = Taylor_values[:5]
    for k in range(3, N + 1):
        x_k = accurate_values[k][0]
        y_k = get_Adams_value(Adams_values, h)
        Adams_values.append([x_k, y_k])
        table.add_row([k, x_k, y_k, abs(y_k - accurate_values[k][1])])
    print("\n\nЭкстраполяционный метод Адамса:")
    print(table)
    return Adams_values[5:]


def Runge_Kutta_method(accurate_values, h, N):
    table = PrettyTable()
    table.field_names = ["k", "x_k", "Приближенное значение y(x_k)", "Абсолютная погрешность"]
    Runge_Kutta_values = [accurate_values[0]]
    for n in range(1, N + 1):
        y_n = Runge_Kutta_values[-1][1]
        x_n = Runge_Kutta_values[-1][0]
        k_1 = h * F.subs({x: x_n, y: y_n})
        k_2 = h * F.subs({x: x_n + h / 2, y: y_n + k_1 / 2})
        k_3 = h * F.subs({x: x_n + h / 2, y: y_n + k_2 / 2})
        k_4 = h * F.subs({x: x_n + h, y: y_n + k_3})
        y_value = y_n + 1 / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4)
        x_value = accurate_values[n][0]
        Runge_Kutta_values.append([x_value, y_value])
        table.add_row([n, x_value, y_value, abs(y_value - accurate_values[n][1])])
    print("\n\nМетод Рунге-Кутта:")
    print(table)
    return Runge_Kutta_values[1:]


def Euler_method(accurate_values, h, N):
    table = PrettyTable()
    table.field_names = ["k", "x_k", "Приближенное значение y(x_k)", "Абсолютная погрешность"]
    Euler_values = [accurate_values[0]]
    for k in range(0, N):
        y_k = Euler_values[-1][1]
        x_k = Euler_values[-1][0]
        y_value = y_k + h * F.subs({x: x_k, y: y_k})
        x_value = accurate_values[k + 1][0]
        Euler_values.append([x_value, y_value])
        table.add_row([k + 1, x_value, y_value, abs(y_value - accurate_values[k + 1][1])])
    print("\n\nМетод Эйлера:")
    print(table)
    return Euler_values[1:]


def Euler_I_method(accurate_values, h, N):
    table = PrettyTable()
    table.field_names = ["k", "x_k", "Приближенное значение y(x_k)", "Абсолютная погрешность"]
    Euler_I_values = [accurate_values[0]]
    for k in range(0, N):
        y_k = Euler_I_values[-1][1]
        x_k = Euler_I_values[-1][0]
        y_k_half = y_k + h / 2 * F.subs({x: x_k, y: y_k})
        y_value = y_k + h * F.subs({x: x_k + h / 2, y: y_k_half})
        x_value = accurate_values[k + 1][0]
        Euler_I_values.append([x_value, y_value])
        table.add_row([k + 1, x_value, y_value, abs(y_value - accurate_values[k + 1][1])])
    print("\n\nМетод Эйлера I:")
    print(table)
    return Euler_I_values[1:]


def Euler_II_method(accurate_values, h, N):
    table = PrettyTable()
    table.field_names = ["k", "x_k", "Приближенное значение y(x_k)", "Абсолютная погрешность"]
    Euler_II_values = [accurate_values[0]]
    for k in range(0, N):
        y_k = Euler_II_values[-1][1]
        x_k = Euler_II_values[-1][0]
        y_k_sup = y_k + h * F.subs({x: x_k, y: y_k})
        x_value = accurate_values[k + 1][0]
        y_value = y_k + h / 2 * (F.subs({x: x_k, y: y_k}) + F.subs({x: x_value, y: y_k_sup}))
        Euler_II_values.append([x_value, y_value])
        table.add_row([k + 1, x_value, y_value, abs(y_value - accurate_values[k + 1][1])])
    print("\n\nМетод Эйлера II:")
    print(table)
    return Euler_II_values[1:]


def task6():
    print("___Численное решение задачи Коши для обычновенного дифференциального уравнения первого порядка___")
    print("\n\nУравнение: y'(x)=" + Y_DERIVATIVE + ", y(0)=" + str(Y_0))
    for i in range(ITER_MAX):
        h = float(input("Введите значение шага h: "))
        N = int(input("Введите количество шагов N: "))

        accurate_solution_table = PrettyTable()
        accurate_solution_table.field_names = ["k", "x_k", "y(x_k)"]
        accurate_values = []
        for k in range(-2, N + 1):
            x_k = round(X_0 + k * h, ROUNDING)
            y_k = Y.subs(x, x_k)
            accurate_solution_table.add_row([k, x_k, y_k])
            accurate_values.append([x_k, y_k])
        print("\n\nТаблица значений точного решения:")
        print(accurate_solution_table)

        Taylor_values = Taylor_series_method(accurate_values, h, N)
        Adams_values = Adams_method(accurate_values[2:], Taylor_values, h, N)  # индексы 0-N
        Runge_Kutta_values = Runge_Kutta_method(accurate_values[2:], h, N)
        Euler_values = Euler_method(accurate_values[2:], h, N)
        Euler_I_values = Euler_I_method(accurate_values[2:], h, N)
        Euler_II_values = Euler_II_method(accurate_values[2:], h, N)

        resp = input('\nХотите продолжить? (Да/Нет): ')
        if resp == 'Нет' or resp == 'нет' or resp == 'No' or resp == 'no':
            sys.exit()
