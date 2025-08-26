// src/gaussder.cpp
#include <cmath>

// Função pulso gaussiano derivado normalizado
inline double gaussder_norm(double t, double m, double sigma) {
    //return (-4.0 / sqrt(2.0 * M_PI)) * ((t - m) / sigma) * std::exp(-(t - m) * (t - m) / (2.0 * sigma * sigma)); //Eq. (3.40)
    return -(std::exp(0.5) / sigma) * (t - m) * std::exp(-std::pow(t - m, 2.0) / (2.0 * sigma * sigma)); //Eq. (3.41)
}
