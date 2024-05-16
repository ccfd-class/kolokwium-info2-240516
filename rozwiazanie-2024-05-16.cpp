#include <cmath>
#include <iostream>

// -----------------------------------------------
// Kod z instrukcji
// -----------------------------------------------

// funkcja liczy wartosc wielomianu interpolacyjnego Lagrange'a
// tablice *x i *y zawieraja wspolrzedne wezlow interpolacji
// n liczba wezlow interpolacji
// xx wartosc dla ktorej liczy sie wielomian
double lagrange(const double *x, const double *y, int n, double xx) {
  int i, j;
  double yint, ylag;

  yint = 0.0;
  for (i = 0; i < n; i++) {
    ylag = 1.0;
    for (j = 0; j < n; j++) {
      if (i == j)
        continue;

      ylag *= (xx - x[j]) / (x[i] - x[j]);
    }

    yint += y[i] * ylag;
  }

  return yint;
}

// oblicza calke metoda simpsona
double simpson(double a, double b, double (*pf)(double), int n) {
  double x = a;
  double h = (b - a) / (2 * n);
  double h2 = h * 2;
  double x1 = a + h;

  double suma = pf(a) + 4. * pf(x1) + pf(b);

  for (int i = 0; i < n - 1; i += 1) {
    x += h2;
    suma += 2. * pf(x) + 4. * pf(x + h);
  }
  return suma * h / 3.;
}

double bisec(double xa, double xb, double (*pf)(double), double eps,
             int *i_iter) {
  int i;
  double fa, fb, xc, fc;

  fa = pf(xa);
  fb = pf(xb);

  if (fa * fb > 0.0) {
    *i_iter = -1;
    return 0;
  }

  for (i = 1; i <= 1000; i++) {
    xc = (xa + xb) / 2.;
    fc = pf(xc);

    if (fa * fc < 0.) {
      xb = xc;
      fb = fc;
    } else {
      xa = xc;
      fa = fc;
    }

    if (fabs(fc) < eps && fabs(xb - xa) < eps)
      break;
  }

  *i_iter = i;
  return xc;
}
// -----------------------------------------------

const double Cd = .2491326;
const double rho = 1025.;
const double A = 80.;
const double vs = 10.;
const double W = 45'000'000.;
const double M = 40'500.;
const double eta = .25;
const int n = 11;
const double delta_x = 200'000;
const double pomiary_vw[] = {3.,  3.5, 3.6, 3.7, 3.,  2.,
                             1.1, .3,  -.4, -1., -1.4};
double pomiary_x[11];

double opor(double x) {
  const double vw = lagrange(pomiary_x, pomiary_vw, n, x);
  const double v_wzgl = vs - vw;
  return .5 * rho * A * Cd * v_wzgl * v_wzgl;
}

double Wa(double d) { return simpson(0., d, opor, 100); }

double bilans_energii(double d) { return Wa(d) - M * W * eta; }

double rozwiazanie() {
  for (int i = 0; i != n; ++i)
    pomiary_x[i] = i * delta_x;
  int n_iter;
  return bisec(0., pomiary_x[n - 1], bilans_energii, 1e-6, &n_iter);
}

int main() {
  const double d = rozwiazanie();
  printf("Obliczony dystans: %.1lf\n", d);
}
