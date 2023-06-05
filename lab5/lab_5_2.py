from math import cos, sin
import numpy as np
from scipy.integrate import quad
from functools import cache
from lab_5_1 import get_polynomial_func, sep, space, get_monom_func, calc_qf


@cache
# Возвращает полиномы Лежандра от 0 до n степени включительно
def calc_gauss_polynoms(n: int) -> list:
    if n == 0:
        return [[1]]
    if n == 1:
        return [[1], [0, 1]]

    polynoms = calc_gauss_polynoms(n - 1)
    # Рекуррентная формула для нахождения полиномов Лежандра
    polynom = np.zeros(n + 1)
    polynom[0] = (1 - n) / n * polynoms[n - 2][0]
    for j in range(n - 2):
        polynom[j + 1] += (2 * n - 1) / n * polynoms[n -
                                                     1][j] + (1 - n) / n * polynoms[n - 2][j + 1]
    for j in range(n - 1, n + 1):
        polynom[j] = (2 * n - 1) / n * polynoms[n - 1][j - 1]
    polynoms.append(list(polynom))
    return polynoms


@cache
def calc_gauss_roots_coeffs(n):
    polynoms = calc_gauss_polynoms(n)
    # Находим корни полиномов Лежандра на [-1, 1]
    roots = list(sorted(np.roots(polynoms[n][::-1])))

    # Находим квадрат предпоследнего полинома
    power_2_polynom_coeffs = [0] * (len(polynoms[n - 1]) * 2 - 1)
    for i in range(len(polynoms[n-1])):
        for j in range(len(polynoms[n-1])):
            power_2_polynom_coeffs[i + j] += polynoms[-2][i] * polynoms[-2][j]
    power_2_polynom = get_polynomial_func(power_2_polynom_coeffs)

    # Находим коэффициенты формулы Гаусса
    coeffs = []
    for i in range(n):
        x = roots[i]
        coeff = 2 * (1 - x ** 2) / (n ** 2 * power_2_polynom(x))
        coeffs.append(coeff)
    return roots, coeffs


# Линейно распространяет корни и коэффициенты на произвольный отрезок
def calc_true_roots_coeffs(a, b, roots, coeffs):
    # Находим коэффициенты формулы Гаусса
    true_roots = []
    true_coeffs = []
    for i in range(len(roots)):
        x = (b - a) / 2 * roots[i] + (b + a) / 2
        coeff = (b - a) / 2 * coeffs[i]
        true_roots.append(x)
        true_coeffs.append(coeff)
    return true_roots, true_coeffs


@cache
def calc_moeller_roots_coeffs(n):
    roots = []
    coeffs = []
    for k in range(1, n + 1):
        roots.append(cos((2 * k - 1) / (2 * n) * np.pi))
        coeffs.append(np.pi / n)
    return roots, coeffs


def gaussian_f_5(x):
    return sin(x)/x


def moeller_f_5(x):
    return cos(x)


if __name__ == "__main__":
    test_func = gaussian_f_5
    moeller_func = moeller_f_5

    print("формула Гаусса")
    sep()
    for n in range(1, 15 + 1):
        print(f"N={n}")
        roots, coeffs = calc_gauss_roots_coeffs(n)
        print("Корни:", " ".join(map(str, roots)))
        print("Коэффициенты:", " ".join(map(str, coeffs)))
        print(
            f"Проверка КФ Гаусса на одночлене x^(2 * {n} - 1) = x({2 * n - 1}):")
        for j in range(2 * n, 2 * n + 1):
            func = get_monom_func(j - 1)
            value = calc_qf(roots, coeffs, func)
            right_value = quad(lambda x: func(x), -1, 1)[0]
            print(f"Проверка на одночлене x^{j - 1}:", value,
                  "\nПравильный ответ:", right_value,
                  "\nПогрешность:", abs(right_value - value), sep="\t")
        sep()

    space()
    print("Формула Гаусса для произвольного отрезка")
    sep()
    print("Введите границы отрезка интегрирования:")
    a, b = map(float, input().split())
    sep()
    for n in range(3, 7):
        print(f"N={n}")
        roots, coeffs = calc_gauss_roots_coeffs(n)
        true_roots, true_coeffs = calc_true_roots_coeffs(a, b, roots, coeffs)
        print("Корни:", " ".join(map(str, true_roots)))
        print("Коэффициенты:", " ".join(map(str, true_coeffs)))
        print(
            f"Проверка КФ Гаусса на одночлене x^(2 * {n} - 1) = x({2 * n - 1}):")
        func = test_func
        # func = get_monom_func(j - 1)
        value = calc_qf(true_roots, true_coeffs, func)
        right_value = quad(lambda x: func(x), a, b)[0]
        print(f"Проверка на тестовой функции для {n} узлов", value,
              "\nПравильный ответ:", right_value,
              "\nПогрешность:", abs(right_value - value), sep="\t")
        sep()

    space()
    print("Формула Мёллера")
    sep()
    print("Введите N1, N2, N3, N4:")
    ns = map(int, input().split())
    sep()
    for n in ns:
        print(f"N={n}")
        roots, coeffs = calc_moeller_roots_coeffs(n)
        print("Корни:", " ".join(map(str, roots)))
        print("Коэффициенты:", " ".join(map(str, coeffs)))
        value = calc_qf(roots, coeffs, moeller_func)
        print("Значение интеграла по КФ Мёллера:", value, sep="\t")
        right_value = quad(lambda x: 1 / (1 - x ** 2) **
                           (1/2) * moeller_func(x), -1, 1)[0]
        print("Значение интеграла по scipy:", right_value, sep="\t")
        print("Погрешность:", abs(right_value - value), sep="\t")
