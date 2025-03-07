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
void simpson_method(double (*)(double), double, double, double, int = 4);

int main()
{
    double a = 1.1;
    double b = 2.631;
    double eps = pow(10, -12);

    trapez_method(f_x, a, b, eps);
    simpson_method(f_x, a, b, eps);
    return 0;
}
double f_x(double x) {
    return (1 + x + pow(x, 2)) / sqrt(pow(x, 3) - 1);
}

void trapez_method(double (*function)(double),double a, double b, double eps, int n) {
    double integral_h = trapez_calc(f_x,a, b, n/2);
    double integral_h2 = trapez_calc(f_x,a, b, n);
    //cout << "n: " << n << endl;
    if (abs(integral_h - integral_h2) <= 3 * eps) {
        cout << setprecision(12) << integral_h2 << endl;
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
    }
    double integral_h = simpson_calc(f_x, a, b, n / 2);
    double integral_h2 = simpson_calc(f_x, a, b, n);
    //cout << "n: " << n << endl;
    if (abs(integral_h - integral_h2) <= 15 * eps) {
        cout << setprecision(12) << integral_h2 << endl;
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