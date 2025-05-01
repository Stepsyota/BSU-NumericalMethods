#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

vector<double> Gauss_solution(vector<vector<double>> A, vector<double> b) {
	int size = A.size();
	vector<int> IOR(size);
	vector<double> x(size);

	for (int i = 0; i < size; ++i) {
		IOR[i] = i;
	}

	for (int k = 0; k < size; ++k) {
		int index_max_el = k;
		double max_value = A[IOR[k]][k];

		for (int i = k + 1; i < size; ++i) {
			if (abs(A[IOR[i]][k]) > abs(max_value)) {
				max_value = A[IOR[i]][k];
				index_max_el = i;
			}
		}
		if (index_max_el != k) {
			swap(IOR[index_max_el], IOR[k]);
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		cout << "1:\n";
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				cout << A[IOR[i]][j] << '\t';
			}
			cout << endl;
		}
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


		if (abs(max_value) < 1e-10) {
			cout << "Matrix has a zero column!\n";
			return vector<double>(size, 0);
		}

		double main_element = max_value;
		for (int j = k; j < size; ++j) {
			A[IOR[k]][j] /= main_element;
		}
		b[IOR[k]] /= main_element;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		cout << "2:\n";
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				cout << A[IOR[i]][j] << '\t';
			}
			cout << endl;
		}
		cout << "b:\n";
		for (int j = 0; j < size; ++j) {
			cout << b[IOR[j]] << '\n';
		}
		cout << endl;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		for (int i = k + 1; i < size; ++i) {
			int str_now = IOR[i];
			double divider = A[str_now][k];

			for (int j = k; j < size; ++j) {
				A[str_now][j] -= divider * A[IOR[k]][j];
			}
			b[str_now] -= divider * b[IOR[k]];
		}
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		cout << "3:\n";
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				cout << A[IOR[i]][j] << '\t';
			}
			cout << endl;
		}
		cout << "b:\n";
		for (int j = 0; j < size; ++j) {
			cout << b[IOR[j]] << '\n';
		}
		cout << endl;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		for (int k = size - 1; k >= 0; --k) {
			x[k] = b[IOR[k]];
			for (int j = k + 1; j < size; ++j) {
				x[k] -= A[IOR[k]][j] * x[j];
			}
		}
	}
	return x;
}

int main()
{
	vector<vector<double>> A =
	{{6, 13, -17},
	{13, 29, -38},
	{-17, -38, 50}};

	vector<double> b = { 2, 4, -5 };
	vector<double> x = Gauss_solution(A, b);
	for (int j = 0; j < x.size(); ++j) {
		cout << x[j] << '\n';
	}
	cout << endl;
}