// f(x) = (1 + x + x^2)/(x^3-1)^(1/2)
// [a;b] = [1.0; 2.631]
// e = 10^(-4), 10^(-5)

#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

double f_x(double);
double trapez_calc(double (*)(double), double, double, int = 1);
void trapez_method(double (*)(double), double, double, double, int = 1);
double simpson_calc(double (*)(double), double, double, int = 4);
void simpson_method(double (*)(double), double, double, double, int = 2);

double f_xy(double, double);
double cubature_simpson_calc(double, double, double, double, int, int);
void cubature_simpson_method(double (*function)(double, double), double, double, double, double, double, int, int);

int main()
{
    double a = 1.1;
    double b = 2.631;
    double eps = pow(10, -5);

    cout << "|" << setw(59) << "The trapezoid method" << setw(39) << "|" << endl;
    cout << "|" << setw(10) << "n1" << setw(10) << "|" << setw(19) << "Integral 1" << setw(10) << "|" << setw(10) << "n2" << setw(10) << "|" << setw(19) << "Integral 2" << setw(10) << "|" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    trapez_method(f_x, a, b, eps);

    cout << "|" << setw(59) << "The Simpson method" << setw(39) << "|" << endl;
    cout << "|" << setw(10) << "n1" << setw(10) << "|" << setw(19) << "Integral 1" << setw(10) << "|" << setw(10) << "n2" << setw(10) << "|" << setw(19) << "Integral 2" << setw(10) << "|" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    simpson_method(f_x, a, b, eps);

    double a1 = 0;
    double b1 = 4.0;
    double c1 = 1.0;
    double d1 = 2.0;
    cout << "|" << setw(59) << "The Simpson Cubature method" << setw(39) << "|" << endl;
    cout << "|" << setw(5) << "N" << setw(5) << "|" << setw(5) << "M" << setw(5) << "|" << setw(19) << "Integral 1" << setw(10) << "|" << setw(5) << "N" << setw(5) << "|" << setw(5) << "M" << setw(5) << "|" << setw(19) << "Integral 2" << setw(10) << "|" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cubature_simpson_method(f_xy, a1, b1, c1, d1, eps, 1, 1);

    return 0;
}

double f_x(double x) {
    return (1 + x + pow(x, 2)) / sqrt(pow(x, 3) - 1);
}

void trapez_method(double (*function)(double),double a, double b, double eps, int n) {
    double integral_h = trapez_calc(f_x,a, b, n/2);
    double integral_h2 = trapez_calc(f_x,a, b, n);
    cout << "|" << setw(10) << n << setw(10) << "|" << setw(19) << integral_h << setw(10) << "|" << setw(10) << 2 * n << setw(10) << "|" << setw(19) << integral_h2 << setw(10) << "|" << endl;
    if (abs(integral_h - integral_h2) <= 3 * eps) {
        cout << "|" << setw(35)  << "The answer: " << setprecision(6) << integral_h2 << " when divided into " << 2 * n <<  " segments" << setw(24) << "|" << endl;
        cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    }
    else {
        trapez_method(f_x, a, b, eps, n * 2);
    }
}

double trapez_calc(double (*function)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double sigma = 0;
    for (int i = 1; i <= n - 1; ++i)
        sigma += function(a + h * i);
    double result = h / 2 * (function(a) + 2 * sigma + function(b));
    return result;
}

void simpson_method(double (*function)(double), double a, double b, double eps, int n) {
    if (n % 2 != 0) {
        cout << "n must be even" << endl;
        exit(1);
    }
    double integral_h = simpson_calc(f_x, a, b, n / 2);
    double integral_h2 = simpson_calc(f_x, a, b, n);
    cout << "|" << setw(10) << n << setw(10) << "|" << setw(19) << integral_h << setw(10) << "|" << setw(10) << 2 * n << setw(10) << "|" << setw(19) << integral_h2 << setw(10) << "|" << endl;
    if (abs(integral_h - integral_h2) <= 15 * eps) {
        cout << "|" << setw(35) << "The answer: " << setprecision(6) << integral_h2 << " when divided into " << 2 * n << " segments" << setw(25) << "|" << endl;
        cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    }
    else {
        simpson_method(f_x, a, b, eps, n * 2);
    }
}

double simpson_calc(double (*function)(double), double a, double b, int n) {
    int m = n/2;
    double h = (b - a) / n;

    double sigma1 = 0;
    for (int i = 1; i <= m; ++i)
        sigma1 += function(a + (2 * i - 1) * h);

    double sigma2 = 0;
    for (int i = 1; i <= m - 1; ++i)
        sigma2 += function(a + 2 * i * h);

    double result = h / 3 * (function(a) + 4 * sigma1 + 2 * sigma2 + function(b));
    return result;
}


double f_xy(double x, double y) {
    return pow(x, 2) / (1 + pow(y, 2));
}

void cubature_simpson_method(double (*function)(double, double), double a, double b, double c, double d, double eps, int N, int M) {
    double integral_h = cubature_simpson_calc(a, b, c, d, N, M);
    double integral_h2 = cubature_simpson_calc(a, b, c, d, 2* N, 2* M);
    cout << "|" << setw(5) << 2 * N << setw(5) << "|" << setw(5) << 2 * M << setw(5) << "|" << setw(19) << integral_h << setw(10) << "|" << setw(5) << 4 * N << setw(5) << "|" << setw(5) << 4 * M << setw(5) << "|" << setw(19) << integral_h2 << setw(10) << "|" << endl;
    if (abs(integral_h - integral_h2) <= 15 * eps) {
        cout << "|" << setw(35) << "The answer: " << setprecision(6) << integral_h2 << " when divided into " << 4* N  << ", " << 4 * M << " segments" << setw(24) << "|" << endl;
        cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    }
    else {
        cubature_simpson_method(f_xy, a, b, c, d, eps, 2*N, 2*M);
    }

}

double cubature_simpson_calc(double a, double b, double c, double d, int N, int M) {
    double hx = (b - a) / (2 * N);
    double hy = (d - c) / (2 * M);

    double sigma = 0;
    for (int i = 0; i <= N - 1; ++i) {
        for (int j = 0; j <= M - 1; ++j) {
            sigma += f_xy(a + (2 * i) * hx, c + (2 * j) * hy); //2i 2j
            sigma += 4 * f_xy(a + (2 * i + 1) * hx, c + (2 * j) * hy); // 2i+1 2j
            sigma += f_xy(a + (2 * i + 2) * hx, c + (2 * j) * hy); // 2i+2 2j
            sigma += 4 * f_xy(a + (2 * i) * hx, c + (2 * j + 1) * hy); //2i 2j+1
            sigma += 16 * f_xy(a + (2 * i + 1) * hx, c + (2 * j + 1) * hy); //2i+1 2j+1
            sigma += 4 * f_xy(a + (2 * i + 2) * hx, c + (2 * j + 1) * hy); //2i+2 2j+1
            sigma += f_xy(a + (2 * i) * hx, c + (2 * j + 2) * hy); //2i 2j+2
            sigma += 4 * f_xy(a + (2 * i + 1) * hx, c + (2 * j + 2) * hy); // 2i+1 2j+2
            sigma += f_xy(a + (2 * i + 2) * hx, c + (2 * j + 2) * hy); // 2i+2 2j+2
        }
    }
    double result = (hx * hy) / 9 * sigma;
    return result;
}