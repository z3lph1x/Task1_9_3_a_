
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

vector<double> Runge_Kutta_4(double y0, double z0, double a, double b, int N) {
    double h = (b - a) / N;
    vector<double> x;
    x.resize(N);
    x[0] = a;
    for (int i = 1; i < N; i++)
        x[i] += i * h;
    vector<double> y;
    y.resize(N);
    vector<double> z;
    z.resize(N);
    y[0] = y0;
    z[0] = z0;
    double k1, k2, k3, k4, l1, l2, l3, l4;
    for (int i = 1; i < N; i++) {
        k1 = h * z[i - 1];
        l1 = h * x[i - 1] * sqrt(y[i - 1]);
        k2 = h * (z[i - 1] + l1 / 2);
        l2 = h * (x[i - 1] + h / 2) * sqrt(y[i - 1] + k1 / 1);
        k3 = h * (z[i - 1] + l2 / 2);
        l3 = h * (x[i - 1] + h / 2) * sqrt(y[i - 1] + k2 / 2);
        k4 = h * (z[i - 1] + l3);
        l4 = h * (x[i - 1] + h / 2) * sqrt(y[i - 1]);
        y[i] = y[i - 1] + 1. / 6 * (k1 + 2 * (k2 + k3) + k4);
        z[i] = z[i - 1] + 1. / 6 * (l1 + 2 * (l2 + l3) + l4);

    }
    return y;
}


int main() {
    int N = 20;
    double e = 0.001;
    double a = 0;
    double b = 1;
    double a1 = 1.5;
    double a2 = 2.5;
    double m = (a1 + a2) / 2;
    vector<double> ans1;
    ans1 = Runge_Kutta_4(0, a1, a, b, N);//alpha1 = a1
    vector<double> ans2 = Runge_Kutta_4(0, a2, a, b, N);//alpha1 = a2
    vector<double> ans3;
    ans3 = Runge_Kutta_4(0, m, a, b, N);//in the middle
    double f1 = ans1[N - 1];
    double f2 = ans2[N - 1];
    double f;
    f = ans3[N - 1];
    while (abs(f - 2) > e) {
        if (f > 2) {
            a2 = m;
            m = (a1 + a2) / 2;
        } else {
            if (b < 2) {
                a1 = m;
                m = (a1 + a2) / 2;
            }
        }
        vector<double> result = Runge_Kutta_4(0, m, a, b, N);
        f = result[N - 1];
    }
    cout << f;
}
