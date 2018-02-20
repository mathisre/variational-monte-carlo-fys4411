#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    double Psi = 1;
    double Phi_k;
    double R_kj;
    double R_kj_2;
    double f;
    for (int k = 0; k < m_system.getNumberOfParticles(); k++){
        Phi_k = 0;
        R_kj_2 = 0;
        f = 1;
        for (int j = k+1; j < m_system->getNumberOfParticles(); j++){
            for (int d = 0; d < m_system.getNumberOfDimensions(); d++){
                R_kj=  particles[j]->getPosition()[d] - particles[k]->getPosition()[d];
                R_kj_2 += R_ij * R_ij;
            }
            R_kj = sqrt(R_kj_2);  //Really only need R_kj, faster way to do this?
            f *= (1 - a / R_kj) * (R_kj > a); //Last part returns 0 or 1 if statement is true (think so)
        }
        for (int d = 0; d < m_system.getNumberOfDimensions(); d++){
            Phi_k +=  (particles[k]->getPosition()[d] * particles[k]->getPosition()[d]); //Add the (1, 1, beta) parameter here
        }
        Phi_k = exp(-getParameters() * Phi_k);
        Psi = Phi_k * f;  //Some way to end loop if f = 0? If test or is that too slow
    }
    return Psi;


    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */
}

double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    double Psi = evaluate(particles);
    // Need to look at the sum of the second derivative.
    // Expression is with the phi's and u's, but evaluate only returns complete wavefunction



    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */
    return 0;
}
