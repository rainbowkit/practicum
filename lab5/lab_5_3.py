from math import sin
from scipy.integrate import quad
from lab_5_2 import calc_gauss_roots_coeffs, calc_true_roots_coeffs
from lab_5_1 import sep, calc_qf


def rho_4(x):
    return x ** 0.5


def calc_true_complex_roots_coeffs(a, b, points, coeffs, m):
    if len(points) != len(coeffs):
        raise ValueError("Количество узлов и коэффициентов должно совпадать")

    eps = (b - a) / m
    true_points = []
    true_coeffs = []
    b_i = a
    for _ in range(m):
        a_i = b_i
        b_i = a_i + eps
        tmp_roots, tmp_coeffs = calc_true_roots_coeffs(
            a_i, b_i, points, coeffs)
        true_points.extend(tmp_roots)
        true_coeffs.extend(tmp_coeffs)
    return true_points, true_coeffs


if __name__ == "__main__":
    test_f = sin
    test_rho = rho_4
    print("Введите a, b - границы отрезка интегрирования")
    a, b = map(float, input().split())

    print("Введите N - количество узлов, m - число отрезков разбиения")
    line = input()
    while line:
        n, m = map(int, line.split())
        roots, coeffs = calc_gauss_roots_coeffs(n)
        true_roots, true_coeffs = calc_true_complex_roots_coeffs(
            a, b, roots, coeffs, m)
        value = calc_qf(true_roots, true_coeffs,
                        lambda x: test_f(x) * test_rho(x))
        right_value = quad(lambda x: test_f(x) * test_rho(x), a, b)[0]
        print(f"Значение интеграла для {n} узлов * {m} разбиений", value,
              "\nПравильный ответ:", right_value,
              "\nПогрешность:", abs(right_value - value), sep="\t")
        sep()
        print("Введите N - количество узлов, m - число отрезков разбиения")
        line = input()
