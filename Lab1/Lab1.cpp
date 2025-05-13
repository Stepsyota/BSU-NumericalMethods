#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

// Прототипы функций
void printMatrix(const vector<vector<double>>& matrix, const vector<double>& b);
vector<double> gaussElimination(vector<vector<double>> A, vector<double> b, double& residualNorm, double& errorNorm);
vector<double> solveLDLT(vector<vector<double>> A, vector<double> b);
void printSolution(const vector<double>& solution, const string& methodName);
double calculateResidualNorm(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b);
vector<double> calculateResidualVector(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b);

int main() {
    setlocale(LC_ALL, "Russian");
    // Пример системы из задачи
    vector<vector<double>> A1 = {
        {2.5, -3.0, 4.6},
        {-3.5, 2.6, 1.5},
        {-6.5, -3.5, 7.3}
    };
    vector<double> b1 = { -1.05, -14.46, -17.73};

    // Пример системы из задачи 21 с λ1=1, λ2=1e3, λ3=1e6
    double lambda1 = 1.0, lambda2 = 1e3, lambda3 = 1e6;
    vector<vector<double>> A2 = {
        {2 * lambda1 + 4 * lambda2, 2 * (lambda1 - lambda2), 2 * (lambda1 - lambda2)},
        {2 * (lambda1 - lambda2), 2 * lambda1 + lambda2 + 3 * lambda3, 2 * lambda1 + lambda2 - 3 * lambda3},
        {2 * (lambda1 - lambda2), 2 * lambda1 + lambda2 - 3 * lambda3, 2 * lambda1 + lambda2 + 3 * lambda3}
    };
    vector<double> b2 = {
        -4 * lambda1 - 2 * lambda2,
        -4 * lambda1 + lambda2 + 9 * lambda3,
        -4 * lambda1 + lambda2 - 9 * lambda3
    };

    // Решение методом Гаусса
    cout << "=============================================\n";
    cout << "Решение системы 1 методом Гаусса:\n";
    cout << "=============================================\n";
    double residualNorm1, errorNorm1;
    vector<double> solutionGauss1 = gaussElimination(A1, b1, residualNorm1, errorNorm1);
    printSolution(solutionGauss1, "Метод Гаусса");

    // Решение системы 2 методом Гаусса
    cout << "\n=============================================\n";
    cout << "Решение системы 2 (задача 21) методом Гаусса:\n";
    cout << "=============================================\n";
    double residualNorm2, errorNorm2;
    vector<double> solutionGauss2 = gaussElimination(A2, b2, residualNorm2, errorNorm2);
    printSolution(solutionGauss2, "Метод Гаусса");

    // Решение системы 2 методом LDLT
    cout << "\n=============================================\n";
    cout << "Решение системы 2 (задача 21) методом LDLT:\n";
    cout << "=============================================\n";
    vector<double> solutionLDLT2 = solveLDLT(A2, b2);
    printSolution(solutionLDLT2, "Метод LDLT");

    return 0;
}

// Функция для печати матрицы и вектора
void printMatrix(const vector<vector<double>>& matrix, const vector<double>& b) {
    int n = matrix.size();
    for (int i = 0; i < n; ++i) {
        cout << "| ";
        for (int j = 0; j < n; ++j) {
            cout << setw(10) << matrix[i][j] << " ";
        }
        cout << " |  | " << setw(10) << b[i] << " |\n";
    }
}

// Функция для печати решения
void printSolution(const vector<double>& solution, const string& methodName) {
    cout << "\n" << methodName << " решение:\n";
    int n = solution.size();
    for (int i = 0; i < n; ++i) {
        cout << "x[" << i + 1 << "] = " << scientific << setprecision(12) << solution[i] << endl;
    }

    cout << fixed << setprecision(6);
}

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
    cout << "\nНорма невязки: " << scientific << residualNorm << endl;
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

    cout << "Относительная погрешность: " << scientific << errorNorm << endl;
    cout << fixed << setprecision(6);
    return x;
}

// Метод LDLT-факторизации
vector<double> solveLDLT(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<double> D(n, 0);

    // Разложение A = LDL^T
    for (int j = 0; j < n; ++j) {
        // Вычисление элементов L и D
        double sum = 0;
        for (int k = 0; k < j; ++k) {
            sum += L[j][k] * L[j][k] * D[k];
        }
        D[j] = A[j][j] - sum;

        for (int i = j + 1; i < n; ++i) {
            sum = 0;
            for (int k = 0; k < j; ++k) {
                sum += L[i][k] * L[j][k] * D[k];
            }
            L[i][j] = (A[i][j] - sum) / D[j];
        }
    }

    // Решение Ly = b (прямая подстановка)
    vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int k = 0; k < i; ++k) {
            sum += L[i][k] * y[k];
        }
        y[i] = (b[i] - sum);
    }

    // Решение Dz = y
    vector<double> z(n);
    for (int i = 0; i < n; ++i) {
        z[i] = y[i] / D[i];
    }

    // Решение L^T x = z (обратная подстановка)
    vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0;
        for (int k = i + 1; k < n; ++k) {
            sum += L[k][i] * x[k];
        }
        x[i] = (z[i] - sum);
    }

    return x;
}

// Вычисление нормы невязки
double calculateResidualNorm(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b) {
    vector<double> residual = calculateResidualVector(A, x, b);
    double norm = 0;
    for (double val : residual) {
        if (fabs(val) > norm) norm = fabs(val);
    }
    cout << "Вектор невязки F:\n";
    for (int i = 0; i < A.size(); ++i) {
        cout << "F[" << i + 1 << "] = " << scientific << residual[i] << endl;
    }
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