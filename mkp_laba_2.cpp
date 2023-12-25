#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>

int testTochnost(double E1, double E2, double tochnost) {
    if (fabs(E1 - E2) <= tochnost) { return 1; }
    else { return 0; }
}

double EkscentrAnomIter(double Ei, double M, double e, double tochnost) {
    double EI = Ei;
    Ei = M + e * sin(EI);

    if (testTochnost(Ei, EI, tochnost) == 1) { return Ei; }
    EkscentrAnomIter(Ei, M, e, tochnost);
}

int main() {

    double R, r_p, r_a, E, M, Mm, tochnost, An, t;
    int stepen;

    std::cout << "Введите: r_p, r_a, R, Mпл *10, Степень,  Точность (0.0001)" << std::endl;
    std::cin >> r_p >> r_a >> R >> Mm >> stepen >> tochnost;

    Mm = Mm * pow(10, stepen);
    double e = (r_a - r_p) / (r_a + r_p + R * 2);

    std::ofstream outTime, outM, outE, outAn;
    outM.precision(200);
    outE.precision(200);
    outTime.precision(200);
    outAn.precision(200);

    outE.open("E.txt", std::ios::app);
    outM.open("M.txt", std::ios::app);
    outTime.open("T.txt", std::ios::app);
    outAn.open("An.txt", std::ios::app);

    std::ofstream outr, outVr, outVn, outV;
    std::cout.precision(200);
    outr.precision(200);
    outVr.precision(200);
    outVn.precision(200);
    outV.precision(200);

    outr.open("r.txt", std::ios::app);
    outVr.open("Vr.txt", std::ios::app);
    outVn.open("Vn.txt", std::ios::app);
    outV.open("V.txt", std::ios::app);

    double G = 6.67428 * pow(10, -20), p = (1.0 / 2.0) * (r_a + r_p + R * 2) * (1 - e * e);
    double Vr, Vn, V, r;
    std::cout << "Радиус нач итер: " << (p / (1 + e * cos(0))) << std::endl << "Скорость нач итер" << pow(G * Mm / p, 0.5) * (1 + e * cos(0));

    if (outAn and outE and outM and outTime) {
        for (t = 0; t <= 25000; t++) {
            M = 2 * M_PI * (t / 25000);
            E = EkscentrAnomIter(M, M, e, tochnost);
            if (E / M_PI <= 1) { An = 2 * atan(pow((1 + e) / (1 - e), 0.5) * tan((E) / 2)); }
            else { An = 2 * M_PI + 2 * atan(pow((1 + e) / (1 - e), 0.5) * tan((E) / 2)); };

            r = (p / (1 + e * cos(An)));
            Vr = pow(G * Mm / p, 0.5) * e * sin(An);
            Vn = pow(G * Mm / p, 0.5) * (1 + e * cos(An));
            V = pow(Vn * Vn + Vr * Vr, 0.5);

            outTime << t / 25000 << std::endl;
            outM << M / M_PI << std::endl;
            outE << E / M_PI << std::endl;
            outAn << An / M_PI << std::endl;

            outr << r << std::endl;
            outVr << Vr << std::endl;
            outVn << Vn << std::endl;
            outV << V << std::endl;


        }
        outTime.close();
        outM.close();
        outE.close();
        outAn.close();

        outr.close();
        outVr.close();
        outVn.close();
        outV.close();
    }
}
