from scipy.integrate import quad
from scipy import linalg
import numpy as np
from math import sin, exp, cos

WEIGHT_MOMENT_ACCURACY = int(1e3)


def rho(x):
    return x**0.5


def f(x):
    return sin(x)


def get_monom_func(n: int):
    def wrapper(x):
        return x ** n
    return wrapper


def get_polynomial_func(coeffs: list):
    def wrapper(x):
        res = 0
        for i in range(len(coeffs)):
            res += coeffs[i] * x ** i
        return res
    return wrapper


def calc_qf_mean_square(a: float, b: float, n: int, func: callable):
    eps = (b-a) / n
    result = 0
    for i in range(n):
        value = eps * func((a + (i + 1/2) * eps))
        result += value
    return result


def calc_weight_moments(a, b, n, rho):
    weight_moments = []
    for i in range(n):
        # TODO изменить на qf_mean_square
        weight_moments.append(quad(lambda x: rho(x) * x ** i, a, b)[0])
    return weight_moments


def calc_iqf_coeffs(a, b, points, rho):  # Интерполяционная квадратурная формула
    n = len(points)
    weight_moments = calc_weight_moments(a, b, n, rho)

    # Находим коэффициенты
    a_matrix = np.zeros((n, n))
    b_vector = np.array(weight_moments)

    # Составляем матрицу коэффициентов
    for i in range(n):
        for j in range(n):
            a_matrix[i][j] = points[j] ** i

    # Решаем систему
    coeffs = linalg.solve(a_matrix, b_vector)
    for i in range(n):
        print(f"Момент веса Mu_{i}:", weight_moments[i],
              f"Коэффициент A_{i + 1}:", coeffs[i],
              f"Узел X_{i + 1}:", points[i], sep="\t")
    return coeffs


def get_hada_polynom(a, b, n, rho):  # Получение ортогонального многочлена степени N=n
    # Получаем весовые моменты
    weight_moments = calc_weight_moments(a, b, 2 * n, rho)

    # Получаем матрицу коэффициентов
    a_matrix = np.zeros((n, n))
    b_vector = -np.array(weight_moments[n:])
    for i in range(n):
        for j in range(n):
            a_matrix[i][j] = weight_moments[i + j]

    # Решаем систему
    coeffs = linalg.solve(a_matrix, b_vector)
    return list(coeffs) + [1]  # Чтобы многочлен был унитарным


# Подсчёт интеграла с помощью КФ НАСТ
def calc_qf_hada_roots_coeffs(a, b, n, rho):
    coeffs = get_hada_polynom(a, b, n, rho)
    
    roots = np.roots(coeffs[::-1])  # Ищем корни многочлена
    roots.sort()

    # # Тест на мономах
    # for i in range(n):
    #     print(quad(lambda x: rho(x) * get_monom_func(i)(x)
    #           * get_polynomial_func(coeffs)(x), a, b))
    # print(roots)

    a_matrix = np.zeros((n, n))
    b_vector = np.array(calc_weight_moments(a, b, n, rho))
    for i in range(n):
        for j in range(n):
            a_matrix[i][j] = roots[j] ** i
    # Решаем систему
    coeffs = linalg.solve(a_matrix, b_vector)
    return roots, coeffs


def calc_qf(points, coeffs, func):
    if len(points) != len(coeffs):
        raise ValueError("Количество узлов и коэффициентов должно совпадать.")
    return sum(coeffs[i] * func(points[i]) for i in range(len(points)))


def sep():
    print("-" * 120)


def space():
    print("\n" * 2)


if __name__ == "__main__":
    test_func = f
    test_rho = rho

    print("Введите a, b:")
    a, b = map(float, input().split())
    print("Введите число узлов N:")
    n = int(input())
    print(f"Введите узлы:")
    tmp_str = input()
    if tmp_str:
        points = list(map(float, tmp_str.split()))
    else:
        points = list(np.linspace(a, b, n))
        print(points)

    sep()
    # Результат scipy
    scipy_value, eps = quad(lambda x: test_func(x) * test_rho(x), a, b)
    print("Результат Scipy: ", scipy_value)

    sep()
    space()
    # ИКФ
    print("ИКФ")
    sep()
    coeffs = calc_iqf_coeffs(a, b, points, test_rho)
    res = calc_qf(points, coeffs, test_func)
    print("Результат ИКФ:", res, "\nПогрешность: ",
          abs(scipy_value - res), sep="\t")

    sep()
    # Проверка ИКФ на одночлене x^(n - 1)
    for j in range(n, n + 1):
        func = get_monom_func(j - 1)
        res = calc_qf(points, coeffs, func)
        right_answer = quad(lambda x: func(x) * test_rho(x), a, b)[0]
        print(f"Проверка ИКФ на одночлене x^{j - 1}:", res,
              "\nПравильный ответ:", right_answer,
              "\nПогрешность:", abs(right_answer - res), sep="\t")

    sep()
    space()
    # КФ НАСТ
    print("КФ НАСТ")
    sep()
    roots, coeffs = calc_qf_hada_roots_coeffs(a, b, n, test_rho)
    print("Корни:", " ".join(map(str, roots)), sep="\t")
    print("Коэффициенты:", " ".join(map(str, coeffs)), sep="\t")
    value = calc_qf(roots, coeffs, test_func)
    print("Значение интеграла по КФ НАСТ:", value,
          "\nПогрешность:", abs(scipy_value - value), sep="\t")

    sep()
    # Проверка КФ НАСТ на одночлене x^(2n - 1)
    print(f"Проверка КФ НАСТ на одночлене x^(2 * {n} - 1) = x({2 * n - 1}):")
    for j in range(2 * n, 2 * n + 1):
        func = get_monom_func(j - 1)
        value = calc_qf(roots, coeffs, func)
        right_answer = quad(lambda x: func(x) * test_rho(x), a, b)[0]
        print(f"Проверка на одночлене x^{j - 1}:", value,
              "\nПравильный ответ:", right_answer,
              "\nПогрешность:", abs(right_answer - value), sep="\t")
    sep()
