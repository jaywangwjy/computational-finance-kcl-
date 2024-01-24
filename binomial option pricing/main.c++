#include <iostream>
#include <cmath>
using namespace std;

double maxi(double x, double y){
    return (x > y) ? x : y;
}

double btPricer(double S0, double K, double r, double T, double vol, int n, bool isAmerican) {
    double value, dt = T / n, nu = r - 0.5 * vol * vol, sqdt = sqrt(dt),
    Delta = sqrt(vol * vol * dt + nu * nu * dt * dt), p = 0.5 + 0.5 * nu * dt / Delta;

    double *St = new double[n + 1];
    double *V = new double[n + 1];
    int i, j;

    // At maturity (i.e. at j == n)
    // X[0] = x0 (i.e. the ln(spot0)) - n * delta
    St[0] = S0 * exp(-n * Delta);
    for (i = 1; i <= n; i++) {
        St[i] = St[i - 1] * exp(2 * Delta);
    }

    for (i = 0; i <= n; i++) {
        V[i] = maxi(St[i] - K, 0);
    }

    // Backward induction
    if (isAmerican) {
        for (j = n - 1; j >= 0; j--) {
            for (i = 0; i <= j; i++) {
                double exercise = maxi(St[i] - K, 0);
                double continuation = exp(-r * dt) * (p * V[i + 1] + (1 - p) * V[i]);
                V[i] = maxi(exercise, continuation);
            }
        }
    } else {
        for (j = n - 1; j >= 0; j--) {
            for (i = 0; i <= j; i++) {
                V[i] = exp(-r * dt) * (p * V[i + 1] + (1 - p) * V[i]);
            }
        }
    }

    value = V[0];

    delete[] St;
    delete[] V;

    return value;
}

int main() {
    double S0 = 100; // Initial stock price
    double K = 105; // Strike price
    double r = 0.05; // Risk-free interest rate
    double T = 1.0; // Time to maturity (in years)
    double vol = 0.2; // Volatility of the underlying stock
    int n = 100; // Number of time steps
    bool isAmerican = true; // Set to true for American option, false for European option

    double optionPrice = btPricer(S0, K, r, T, vol, n, isAmerican);

    cout << "Option Price: " << optionPrice << endl;

    return 0;
}

