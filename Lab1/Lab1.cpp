#include <iostream>

using namespace std;

const int n = 3;
int main()
{
    int A[3][3]{ {6, 13, -17}, {13, 29, -38}, {-17, -38, 50} };
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            cout << A[i][j] << '\t';
        }
        cout << endl;
    }
    int b[3][1]{ {2}, {4}, {-5} };
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 1; ++j) {
            cout << b[i][j] << '\t';
        }
        cout << endl;
    }

    //int IOR[3];
    //for (int i = 0; i < n; ++i) {
    //    IOR[i] = i;
    //}
    //for (int i = 0; i < n; ++i) {
    //    cout << IOR[i] << '\t';
    //}
    //cout << endl;

    //for (int i = 0; i < n; ++i) {
    //    int AKK = 0;
    //    int num_str_AKK = i;
    //    for (int j = 0; j < n; ++j) {
    //        int index = IOR[j];
    //        if (abs(A[index][i]) > AKK) {
    //            AKK = abs(A[index][i]);
    //            num_str_AKK = i;
    //            cout << "AKK " << AKK << endl;
    //            cout << "i " << num_str_AKK << "\tj " << j << endl;
    //        }
    //    }

    //    if (AKK == 0) {
    //        cout << "Error. Matrix str 0" << endl;
    //        exit(1);
    //    }

    //    swap(IOR[i], IOR[num_str_AKK]);
    // }
    gauss(A, b);

    return 0;
}

void gauss(int a[][n], int b[][1]) {
    int IOR[3];
    for (int i = 0; i < n; ++i) {
        IOR[i] = i;
    }

    // Прямой ход
    for (int k = 0; k < n; ++k) {
        // Поиск главного элемента
        double AKK = 0.0;
        int p = k;
        for (int i = k; i < n; ++i) {
            int index = IOR[i];
            if (fabs(a[index][k]) > AKK) {
                AKK = fabs(a[index][k]);
                p = i;
            }
        }

        // Проверка на вырождение матрицы
        if (AKK == 0.0) {
            cout << "Matrix is singular!" << endl;
            exit(1);
        }

        // Перестановка строк
        swap(IOR[k], IOR[p]);

        // Выбор ведущего элемента
        int lead = IOR[k];
        double AMAIN = a[lead][k];

        // Исключение переменной x_k
        for (int j = k; j < n; ++j) {
            a[lead][j] /= AMAIN;
        }
        b[lead] /= AMAIN;

        for (int i = k + 1; i < n; ++i) {
            int row = IOR[i];
            double multiplier = a[row][k];
            for (int j = k; j < n; ++j) {
                a[row][j] -= multiplier * a[lead][j];
            }
            b[row] -= multiplier * b[lead];
        }
    }
}
