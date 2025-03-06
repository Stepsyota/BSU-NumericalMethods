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

int main()
{
    double a = 1.1;
    double b = 2.631;
    double eps = 12;
   
    trapez_method(f_x, a, b, eps);
    return 0;
}
double f_x(double x) {
    return (1 + x + pow(x, 2)) / sqrt(pow(x, 3) - 1);
}

void trapez_method(double (*function)(double),double a, double b, double eps, int n) {
    double integral_h = trapez_calc(f_x,a, b, n/2);
    double integral_h2 = trapez_calc(f_x,a, b, n);
    if (abs(integral_h - integral_h2) <= 3 * pow(10, -eps)) {
        cout << setprecision(eps) << integral_h2 << endl;
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