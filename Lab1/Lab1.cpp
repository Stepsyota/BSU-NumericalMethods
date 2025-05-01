#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
vector<double> Gauss_solution(vector<vector<double>>, vector<double>);

int main()
{
	vector<vector<double>> A =
	{{6, 13, -17},
	{13, 29, -38},
	{-17, -38, 50}};

	vector<double> b = { 2, 4, -5 };
	vector<double> x = Gauss_solution(A, b);
	for (int j = 0; j < x.size(); ++j) {
		cout << "x" << j +1 << " = " << x[j] << '\n';
	}
	cout << endl;
}

vector<double> Gauss_solution(vector<vector<double>> A, vector<double> b) {

	if (A.size() != b.size()) {
		cout << "ERROR! Invalid sizes\n";
		exit(-1);
	}

	int size = A.size();
	vector<int> IOR(size);
	vector<double> x(size);

	// Array of row order
	for (int i = 0; i < size; ++i) {
		IOR[i] = i;
	}

	//To a triangular form
	for (int k = 0; k < size; ++k) {
		int index_max_el = k;
		double max_value = A[IOR[k]][k];

		//Find max element in a column
		for (int i = k + 1; i < size; ++i) {
			if (abs(A[IOR[i]][k]) > abs(max_value)) {
				max_value = A[IOR[i]][k];
				index_max_el = i;
			}
		}
		// Changing the order of the rows so that the maximum element is the highest.
		if (index_max_el != k) {
			swap(IOR[index_max_el], IOR[k]);
		}

		// Checking for a non-zero column
		if (abs(max_value) < 1e-10) {
			cout << "Matrix has a zero column!\n";
			return vector<double>(size, 0);
		}

		// Normalize the row
		double main_element = max_value;
		for (int j = k; j < size; ++j) {
			A[IOR[k]][j] /= main_element;
		}
		b[IOR[k]] /= main_element;

		// Counting the remaining elements
		for (int i = k + 1; i < size; ++i) {
			int str_now = IOR[i];
			double divider = A[str_now][k];

			for (int j = k; j < size; ++j) {
				A[str_now][j] -= divider * A[IOR[k]][j];
			}
			b[str_now] -= divider * b[IOR[k]];
		}

		// Counting the values of x
		for (int k = size - 1; k >= 0; --k) {
			x[k] = b[IOR[k]];
			for (int j = k + 1; j < size; ++j) {
				x[k] -= A[IOR[k]][j] * x[j];
			}
		}
	}
	return x;
}