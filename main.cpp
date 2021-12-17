
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
        x[i] = x[i-1] + h;
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
    int N = 1000;
    double e = 0.001;
    double a = 0;
    double b = 1;
    double h = (b-a)/N;
    vector<double> alpha;
    alpha.resize(N);
    alpha[0] = 5;
    double ans1 = Runge_Kutta_4(0, alpha[0]+ h, a, b, N)[N-1];//alpha1 + h
    double ans2 = Runge_Kutta_4(0, alpha[0], a, b, N)[N-1];//alpha
    double df = (ans1 - ans2)/h;
    alpha[1]=alpha[0] - (ans2 - 2)/df;
    for(int i = 2;i < N;i++){
        if(abs(alpha[i-1]-alpha[i-2])>e) {
            ans1 = Runge_Kutta_4(0, alpha[i - 1] + h, a, b, N)[N - 1];//alpha1 + h
            ans2 = Runge_Kutta_4(0, alpha[i - 1], a, b, N)[N - 1];//alpha
            df = (ans1 - ans2) / h;
            alpha[i] = alpha[i - 1] - (ans2 - 2) / df;
        }
        else {
            cout<< alpha[i-1]<<"<- alpha; "<<ans2<<" <- definition of function at the last iteration";
            return 0;
        }
    }
}
