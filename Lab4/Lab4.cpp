#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// Прототипы функций
vector<double> gaussElimination(vector<vector<double>> matrix, vector<double> constants);
void performLeastSquaresApproximation(const vector<double>& H, const vector<double>& mu);
void displayResults(const vector<double>& coefficients, const vector<double>& H,
    const vector<double>& mu, double residualVariance);

int main() {
    setlocale(LC_ALL, "Russian");
    // Исходные данные из задания
    vector<double> H_data = { 0.164, 0.328, 0.656, 0.984, 1.312, 1.640 };
    vector<double> mu_data = { 0.448, 0.432, 0.421, 0.417, 0.414, 0.412 };

    // Выполнение аппроксимации методом наименьших квадратов
    performLeastSquaresApproximation(H_data, mu_data);

    return 0;
}

// Функция для решения системы линейных уравнений методом Гаусса
vector<double> gaussElimination(vector<vector<double>> A, vector<double> b) {
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
    return x;
}

// Функция для выполнения аппроксимации методом наименьших квадратов
void performLeastSquaresApproximation(const vector<double>& H, const vector<double>& mu) {
    int dataPointsCount = H.size();

    // Инициализация матрицы системы и вектора правых частей
    vector<vector<double>> normalMatrix(2, vector<double>(2, 0.0));
    vector<double> rightSideVector(2, 0.0);

    // Вычисление сумм для нормальных уравнений
    double sumInvH2 = 0, sumInvH = 0;
    double sumMuInvH = 0, sumMu = 0;

    for (int i = 0; i < dataPointsCount; i++) {
        double invH = 1.0 / H[i];

        sumInvH2 += invH * invH;
        sumInvH += invH;
        sumMuInvH += mu[i] * invH;
        sumMu += mu[i];
    }

    // Заполнение матрицы нормальных уравнений
    normalMatrix[0][0] = sumInvH2;  normalMatrix[0][1] = sumInvH;
    normalMatrix[1][0] = sumInvH;   normalMatrix[1][1] = dataPointsCount;

    // Заполнение вектора правых частей
    rightSideVector[0] = sumMuInvH;
    rightSideVector[1] = sumMu;

    // Решение системы уравнений
    vector<double> coefficients = gaussElimination(normalMatrix, rightSideVector);

    // Вычисление остаточной дисперсии
    double residualSum = 0.0;
    for (int i = 0; i < dataPointsCount; i++) {
        double predictedMu = coefficients[0] / H[i] + coefficients[1];
        residualSum += pow(mu[i] - predictedMu, 2);
    }
    double residualVariance = residualSum / (dataPointsCount - 2); // Число степеней свободы

    // Вывод результатов
    displayResults(coefficients, H, mu, residualVariance);
}

// Функция для вывода результатов
void displayResults(const vector<double>& coeffs, const vector<double>& H,
    const vector<double>& mu, double residualVar) {
    cout << "========================================================" << endl;
    cout << " РЕЗУЛЬТАТЫ АППРОКСИМАЦИИ ЗАВИСИМОСТИ Mu = a/H + b" << endl;
    cout << "========================================================" << endl;

    cout << "\n Полученные коэффициенты модели:" << endl;
    cout << fixed << setprecision(6);
    cout << " a (коэффициент при 1/H) = " << setw(12) << coeffs[0] << endl;
    cout << " b (свободный член) = " << setw(19) << coeffs[1] << endl;

    cout << "========================================================" << endl;
    cout << "\n Статистические показатели качества аппроксимации:" << endl;
    cout << " Остаточная дисперсия: " << setw(26) << residualVar << endl;
    cout << " Среднеквадратичное отклонение: " << setw(16) << sqrt(residualVar) << endl;

    cout << "========================================================" << endl;
    cout << "\n ТАБЛИЦА СРАВНЕНИЯ ИСХОДНЫХ И РАСЧЕТНЫХ ЗНАЧЕНИЙ\n" << endl;
    cout << setw(10) << "H" << setw(15) << "Исходное mu"
        << setw(15) << "Расчетное mu" << setw(15) << "Отклонение" << endl;

    for (size_t i = 0; i < H.size(); i++) {
        double predictedMu = coeffs[0] / H[i] + coeffs[1];
        cout << fixed << setprecision(3);
        cout << setw(10) << H[i];
        cout << fixed << setprecision(6);
        cout << setw(15) << mu[i] << setw(15) << predictedMu
            << setw(15) << (mu[i] - predictedMu) << endl;
    }
    cout << "========================================================" << endl;
}