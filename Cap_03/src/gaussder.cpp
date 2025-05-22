// src/gaussder.cpp
#include <cmath>

// Função gaussiana derivada normalizada
inline double gaussder_norm(double t, double m, double sigma) {
    //return (-1.0 / sqrt(2.0 * M_PI)) * ((t - m) / (sigma * sigma * sigma)) * exp(-(t - m) * (t - m) / (2.0 * sigma * sigma)); //Eq. (3.40)
    return -(std::exp(0.5) / sigma) * (t - m) * std::exp(-std::pow(t - m, 2) / (2.0 * sigma * sigma)); //Eq. (3.41)
}
