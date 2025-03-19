#include <iostream>

using namespace std;

int main()
{
    const int A[3][3]{ {6, 13, -17}, {13, 29, -38}, {-17, -38, 50} };
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            cout << A[i][j] << '\t';
        }
        cout << endl;
    }
    const int b[3][1]{ {2}, {4}, {-5} };
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 1; ++j) {
            cout << b[i][j] << '\t';
        }
        cout << endl;
    }

    int IOR[3];
    for (int i = 0; i < 3; ++i) {
        IOR[i] = i;
    }
    for (int i = 0; i < 3; ++i) {
        cout << IOR[i] << '\t';
    }
    cout << endl;

    for (int i = 0; i < 3; ++i) {
        double AKK = 0;
        int p = i;
        for (int j = i; j < 3; ++j) {
            int index = IOR[j];
            if (abs(A[index][i]) > AKK) {
                AKK = abs(A[index][i]);
                p = i;
                cout << "AKK " << AKK << endl;
                cout << "i " << p << "\tj " << j << endl;
            }
        }
    }
    return 0;
}
