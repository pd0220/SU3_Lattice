#pragma once
#pragma warning(push)
#pragma warning(disable : 4127)

// used headers and/or libraries
// custom
#include "auxiliary.hh"
#include "observables.hh"

using namespace std;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// quark contribution to the determinant at a single site
complex<double> Mq_x(double const &prf, complex<double> const &Px)
{
    // return M_x^(q)
    return 1. + prf * Px + sq(prf) * conj(Px) + cb(prf);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// antiquark contribution to the determinant at a single site
complex<double> MqBar_x(double const &prf, complex<double> const &Px)
{
    // return M_x^(qBar)
    return 1. + prf * conj(Px) + sq(prf) * Px + cb(prf);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// change in phase quenched fermionic action ~ detM -> |detM|
double DeltaFermionAction_PQ(int const &Nf, double const &q_prf, double const &qBar_prf, Eigen::Matrix3cd const &X, vector<Eigen::Matrix3cd> const &links, vector<int> const &coordinates, vector<int> const &dims)
{
    // no change if link is in spatial direction
    if (coordinates[4] != 3)
        return 0.;
    else
    {
        // compute Polyakov loop at given spatial site
        // without taking the trace
        Eigen::Matrix3cd Px = PolyakovLoop_noTrace(links, dims, coordinates);
        // with new link
        Eigen::Matrix3cd Px_new = X * Px;

        // compute differences of the phase quenched action
        double deltaFermionAction_q = log(abs(Mq_x(q_prf, Px_new.trace()))) - log(abs(Mq_x(q_prf, Px.trace())));
        double deltaFermionAction_qBar = log(abs(MqBar_x(qBar_prf, Px_new.trace()))) - log(abs(MqBar_x(qBar_prf, Px.trace())));

        // return
        return -(double)Nf * 2. * (deltaFermionAction_q + deltaFermionAction_qBar);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// fermion determinant (Nf = 1)
complex<double> FermionDeterminant(double const &q_prf, double const &qBar_prf, vector<Eigen::Matrix3cd> const &links, vector<int> const &dims)
{
    // lattices size
    int Ns = dims[0];

    // initialize product
    complex<double> prod = 1.;
    // loop over spatial directions
    for (int x = 0; x < Ns; x++)
        for (int y = 0; y < Ns; y++)
            for (int z = 0; z < Ns; z++)
            {
                // compute Polyakov loop
                vector<int> spatialSite = {x, y, z};
                complex<double> Px = PolyakovLoop(links, dims, spatialSite);
                // contributions to the determinant
                complex<double> Mq = Mq_x(q_prf, Px);
                complex<double> MqBar = MqBar_x(qBar_prf, Px);

                // update product
                prod *= sq(Mq * MqBar);
            }

    // return HDQCD determinant
    return prod;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// density
complex<double> Density(int const &Nf, double const &q_prf, double const &qBar_prf, vector<Eigen::Matrix3cd> const &links, vector<int> const &dims)
{
    // lattice size
    int Ns = dims[0];

    // initialize sum
    complex<double> sum = 0.;
    // loop over spatial directions
    for (int x = 0; x < Ns; x++)
        for (int y = 0; y < Ns; y++)
            for (int z = 0; z < Ns; z++)
            {
                // compute density contribution for a given site
                vector<int> spatialSite = {x, y, z};
                complex<double> Px = PolyakovLoop(links, dims, spatialSite);
                complex<double> Mq = Mq_x(q_prf, Px);
                complex<double> MqBar = MqBar_x(qBar_prf, Px);

                sum += (q_prf * Px / 3. + 2. / 3. * sq(q_prf) * conj(Px) + cb(q_prf)) / Mq - (qBar_prf * conj(Px) / 3. + 2. / 3. * sq(qBar_prf) * Px + cb(qBar_prf)) / MqBar;
            }

    // return density
    return (double)Nf * 6. / (double)cb(Ns) * sum;
}

#pragma warning(pop)