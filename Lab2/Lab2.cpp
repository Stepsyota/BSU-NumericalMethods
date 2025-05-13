#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept>


using namespace std;

// Прототипы функций из первой лабораторной работы
vector<double> gaussElimination(vector<vector<double>> A, vector<double> b, double& residualNorm, double& errorNorm);
double calculateResidualNorm(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b);
vector<double> calculateResidualVector(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b);

// Прототипы новых функций
void printNewtonHeader();
void printNewtonIteration(int iter, double residualNorm, double solutionChangeNorm);
vector<double> solveNewtonAnalytic(
    vector<double>(*systemFunc)(const vector<double>&),
    vector<vector<double>>(*jacobianFunc)(const vector<double>&),
    const vector<double>& initialGuess,
    double epsilon1,
    double epsilon2,
    int maxIterations);
vector<double> solveNewtonNumerical(
    vector<double>(*systemFunc)(const vector<double>&),
    const vector<double>& initialGuess,
    double epsilon1,
    double epsilon2,
    int maxIterations,
    double M);
double vectorNorm(const vector<double>& vec);
double vectorDiffNorm(const vector<double>& vec1, const vector<double>& vec2);

// Пример системы уравнений
vector<double> exampleSystem1(const vector<double>& x) {
    vector<double> F(2);
    F[0] = 1.5 * pow(x[0],3) - pow(x[1], 2) - 1;
    F[1] = x[0] * pow(x[1], 3) - x[1] - 4;
    return F;
}

// Матрица Якоби
vector<vector<double>> exampleJacobian1(const vector<double>& x) {
    vector<vector<double>> J(2, vector<double>(2));
    J[0][0] = 4.5 * pow(x[0], 2);
    J[0][1] = -2 * x[1];
    J[1][0] = pow(x[1], 3);
    J[1][1] = 3 * x[0] * pow(x[1], 2) - 1;
    return J;
}

int main() {
    setlocale(LC_ALL, "Russian");

    // Параметры метода
    const double epsilon1 = 1e-9;
    const double epsilon2 = 1e-9;
    const int maxIterations = 100;

    // Начальное приближение для задачи 1
    vector<double> initialGuess = { 1.0, 1.0 };

    // Решение с аналитическим Якобианом
    cout << "=============================================\n";
    cout << "Решение системы методом Ньютона (аналитический Якобиан):\n";
    cout << "=============================================\n";
    vector<double> solutionAnalytic = solveNewtonAnalytic(
        exampleSystem1, exampleJacobian1, initialGuess, epsilon1, epsilon2, maxIterations);

    cout << "\nРешение:\n";
    for (size_t i = 0; i < solutionAnalytic.size(); ++i) {
        cout << "x[" << i + 1 << "] = " << scientific << setprecision(12) << solutionAnalytic[i] << endl;
    }

    // Решение с численным Якобианом (разные значения параметра M)
    vector<double> M_values = { 0.01, 0.05, 0.1 };
    for (double M : M_values) {
        cout << "\n=============================================\n";
        cout << "Решение системы методом Ньютона (численный Якобиан, M = " << M << "):\n";
        cout << "=============================================\n";

        vector<double> solutionNumerical = solveNewtonNumerical(
            exampleSystem1, initialGuess, epsilon1, epsilon2, maxIterations, M);

        cout << "\nРешение:\n";
        for (size_t i = 0; i < solutionNumerical.size(); ++i) {
            cout << "x[" << i + 1 << "] = " << scientific << setprecision(12) << solutionNumerical[i] << endl;
        }

        // Сравнение с аналитическим решением
        vector<double> diff(solutionNumerical.size());
        for (size_t i = 0; i < diff.size(); ++i) {
            diff[i] = solutionNumerical[i] - solutionAnalytic[i];
        }
        cout << "Разница с аналитическим решением: " << scientific << setprecision(4)
            << vectorNorm(diff) << endl;
    }

    return 0;
}

// Функция для вычисления нормы вектора (максимального по модулю элемента)
double vectorNorm(const vector<double>& vec) {
    double norm = 0.0;
    for (double val : vec) {
        if (abs(val) > norm) {
            norm = abs(val);
        }
    }
    return norm;
}

// Функция для вычисления нормы разности векторов
double vectorDiffNorm(const vector<double>& vec1, const vector<double>& vec2) {
    if (vec1.size() != vec2.size()) {
        throw invalid_argument("Векторы должны быть одного размера");
    }

    double norm = 0.0;
    for (size_t i = 0; i < vec1.size(); ++i) {
        double diff = abs(vec1[i] - vec2[i]);
        // Для значений >= 1 используем относительную разность, иначе абсолютную
        if (abs(vec1[i]) >= 1.0) {
            diff /= abs(vec1[i]);
        }
        if (diff > norm) {
            norm = diff;
        }
    }
    return norm;
}

// Печать заголовка таблицы итераций
void printNewtonHeader() {
    cout << "Iter | Residual norm | Solution change norm" << endl;
    cout << "-------------------------------------------" << endl;
}

// Печать информации о текущей итерации
void printNewtonIteration(int iter, double residualNorm, double solutionChangeNorm) {
    cout << setw(4) << iter << " | " << scientific << setprecision(4)
        << residualNorm << " | " << solutionChangeNorm << endl;
}

// Метод Ньютона с аналитическим Якобианом
vector<double> solveNewtonAnalytic(
    vector<double>(*systemFunc)(const vector<double>&),
    vector<vector<double>>(*jacobianFunc)(const vector<double>&),
    const vector<double>& initialGuess,
    double epsilon1,
    double epsilon2,
    int maxIterations) {

    vector<double> x = initialGuess;
    vector<double> x_prev = x;
    int iter = 0;

    printNewtonHeader();

    while (true) {
        // Вычисление функции системы в текущей точке
        vector<double> F = systemFunc(x);
        double residualNorm = vectorNorm(F);

        // Вычисление матрицы Якоби
        vector<vector<double>> J = jacobianFunc(x);

        // Решение системы J*dx = -F методом Гаусса
        double dummy1, dummy2; // нам не нужны нормы невязки для этой вспомогательной системы
        vector<double> dx = gaussElimination(J, F, dummy1, dummy2);
        for (double& val : dx) val = -val; // так как решаем J*dx = -F

        // Обновление решения
        x_prev = x;
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] += dx[i];
        }

        // Вычисление нормы изменения решения
        double solutionChangeNorm = vectorDiffNorm(x, x_prev);

        // Вывод информации о текущей итерации
        printNewtonIteration(iter, residualNorm, solutionChangeNorm);

        // Проверка критериев остановки
        if (residualNorm <= epsilon1 && solutionChangeNorm <= epsilon2) {
            cout << "Решение найдено с требуемой точностью!" << endl;
            break;
        }

        if (iter >= maxIterations) {
            cout << "Достигнуто максимальное число итераций!" << endl;
            break;
        }

        iter++;
    }

    return x;
}

// Метод Ньютона с численным Якобианом
vector<double> solveNewtonNumerical(
    vector<double>(*systemFunc)(const vector<double>&),
    const vector<double>& initialGuess,
    double epsilon1,
    double epsilon2,
    int maxIterations,
    double M) {

    vector<double> x = initialGuess;
    vector<double> x_prev = x;
    int iter = 0;

    printNewtonHeader();

    while (true) {
        // Вычисление функции системы в текущей точке
        vector<double> F = systemFunc(x);
        double residualNorm = vectorNorm(F);

        // Численное вычисление матрицы Якоби
        vector<vector<double>> J(F.size(), vector<double>(x.size()));
        vector<double> fx = systemFunc(x);

        for (size_t j = 0; j < x.size(); ++j) {
            vector<double> x_perturbed = x;
            double delta = M * x[j];
            if (abs(delta) < 1e-12) delta = M; // Защита от нулевого приращения
            x_perturbed[j] += delta;

            vector<double> fx_perturbed = systemFunc(x_perturbed);

            for (size_t i = 0; i < F.size(); ++i) {
                J[i][j] = (fx_perturbed[i] - fx[i]) / delta;
            }
        }

        // Решение системы J*dx = -F методом Гаусса
        double dummy1, dummy2; // нам не нужны нормы невязки для этой вспомогательной системы
        vector<double> dx = gaussElimination(J, F, dummy1, dummy2);
        for (double& val : dx) val = -val; // так как решаем J*dx = -F

        // Обновление решения
        x_prev = x;
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] += dx[i];
        }

        // Вычисление нормы изменения решения
        double solutionChangeNorm = vectorDiffNorm(x, x_prev);

        // Вывод информации о текущей итерации
        printNewtonIteration(iter, residualNorm, solutionChangeNorm);

        // Проверка критериев остановки
        if (residualNorm <= epsilon1 && solutionChangeNorm <= epsilon2) {
            cout << "Решение найдено с требуемой точностью!" << endl;
            break;
        }

        if (iter >= maxIterations) {
            cout << "Достигнуто максимальное число итераций!" << endl;
            break;
        }

        iter++;
    }

    return x;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Метод Гаусса с выбором главного элемента
vector<double> gaussElimination(vector<vector<double>> A, vector<double> b, double& residualNorm, double& errorNorm) {
    int n = A.size();
    vector<int> rowOrder(n);
    for (int i = 0; i < n; ++i) rowOrder[i] = i;

    // Сохраняем исходные матрицы для вычисления невязки
    vector<vector<double>> originalA = A;
    vector<double> originalB = b;

    // Прямой ход
    for (int k = 0; k < n; ++k) {
        // Поиск главного элемента в столбце k
        int maxRow = k;
        double maxVal = fabs(A[rowOrder[k]][k]);
        for (int i = k + 1; i < n; ++i) {
            if (fabs(A[rowOrder[i]][k]) > maxVal) {
                maxVal = fabs(A[rowOrder[i]][k]);
                maxRow = i;
            }
        }

        // Перестановка строк
        swap(rowOrder[k], rowOrder[maxRow]);

        // Проверка на вырожденность
        if (fabs(A[rowOrder[k]][k]) < 1e-15) {
            cerr << "Матрица вырождена!" << endl;
            exit(1);
        }

        // Нормализация строки k
        double pivot = A[rowOrder[k]][k];
        for (int j = k; j < n; ++j) {
            A[rowOrder[k]][j] /= pivot;
        }
        b[rowOrder[k]] /= pivot;

        // Исключение переменной x_k из остальных уравнений
        for (int i = k + 1; i < n; ++i) {
            double factor = A[rowOrder[i]][k];
            for (int j = k; j < n; ++j) {
                A[rowOrder[i]][j] -= factor * A[rowOrder[k]][j];
            }
            b[rowOrder[i]] -= factor * b[rowOrder[k]];
        }
    }

    // Обратный ход
    vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[rowOrder[i]];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= A[rowOrder[i]][j] * x[j];
        }
    }

    // Вычисление невязки
    residualNorm = calculateResidualNorm(originalA, x, originalB);
    //cout << "\nНорма невязки: " << scientific << residualNorm << endl;
    // Оценка погрешности через вспомогательную систему
    vector<double> Ax(n);
    for (int i = 0; i < n; ++i) {
        Ax[i] = 0;
        for (int j = 0; j < n; ++j) {
            Ax[i] += originalA[i][j] * x[j];
        }
    }

    vector<double> xTilde = x; // начальное приближение
    for (int iter = 0; iter < 5; ++iter) { // несколько итераций уточнения
        vector<double> residual = calculateResidualVector(originalA, xTilde, Ax);
        for (int i = 0; i < n; ++i) {
            xTilde[i] -= residual[i] / originalA[i][i]; // простое итерационное уточнение
        }
    }

    // Вычисление относительной погрешности
    errorNorm = 0;
    double maxX = 0;
    for (int i = 0; i < n; ++i) {
        double diff = fabs(x[i] - xTilde[i]);
        if (diff > errorNorm) errorNorm = diff;
        if (fabs(x[i]) > maxX) maxX = fabs(x[i]);
    }
    errorNorm /= maxX;

    //cout << "Относительная погрешность: " << scientific << errorNorm << endl;
    //cout << fixed << setprecision(6);
    return x;
}

// Вычисление нормы невязки
double calculateResidualNorm(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b) {
    vector<double> residual = calculateResidualVector(A, x, b);
    double norm = 0;
    for (double val : residual) {
        if (fabs(val) > norm) norm = fabs(val);
    }
    //cout << "Вектор невязки F:\n";
    //for (int i = 0; i < A.size(); ++i) {
    //    cout << "F[" << i + 1 << "] = " << scientific << residual[i] << endl;
    //}
    return norm;
}

// Вычисление вектора невязки
vector<double> calculateResidualVector(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b) {
    int n = A.size();
    vector<double> residual(n);
    for (int i = 0; i < n; ++i) {
        residual[i] = -b[i];
        for (int j = 0; j < n; ++j) {
            residual[i] += A[i][j] * x[j];
        }
    }
    return residual;
}