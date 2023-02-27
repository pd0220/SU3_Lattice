#pragma once
#pragma warning(push)
#pragma warning(disable : 4127)

// used headers and/or libraries
// custom
#include "auxiliary.hh"

using namespace std;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// plaquette at given site with given mu and nu Lorentz indices
Eigen::Matrix3cd Plaquette(vector<Eigen::Matrix3cd> const &links, vector<int> const &coordinates, vector<int> const &dims, int const &nu)
{
    // given mu direction
    int mu = coordinates[4];

    // used coordinate combinations in the plaquettes
    //
    // U_mu(n)
    vector<int> n_at_mu = coordinates;
    // U_nu(n + mu)
    vector<int> n_plus_mu_at_nu = coordinates;
    n_plus_mu_at_nu[mu] += 1;
    n_plus_mu_at_nu[4] = nu;
    // U_mu(n + nu)
    vector<int> n_plus_nu_at_mu = coordinates;
    n_plus_nu_at_mu[nu] += 1;
    // U_nu(n)
    vector<int> n_at_nu = coordinates;
    n_at_nu[4] = nu;

    // plaquette
    return links[GetIndex(n_at_mu, dims)] * links[GetIndex(n_plus_mu_at_nu, dims)] * links[GetIndex(n_plus_nu_at_mu, dims)].adjoint() * links[GetIndex(n_at_nu, dims)].adjoint();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// average plaquette
double AveragePlaquette(vector<Eigen::Matrix3cd> const &links, vector<int> const &dims)
{
    // lattices sizes
    int Ns = dims[0];
    int Nt = dims[3];
    int NSite = cb(Ns) * Nt;

    // sum of plaquettes
    double sum = 0.;
    // loop over spacetime coordinates & directions
    for (int x = 0; x < Ns; x++)
        for (int y = 0; y < Ns; y++)
            for (int z = 0; z < Ns; z++)
                for (int t = 0; t < Nt; t++)
                    for (int mu = 0; mu < 4; mu++)
                        for (int nu = mu + 1; nu < 4; nu++)
                        {
                            vector<int> coordinates = {x, y, z, t, mu};
                            sum += Plaquette(links, coordinates, dims, nu).trace().real();
                        }
    // returing normalized average plaquette
    return sum / 18. / (double)NSite;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// average spatial plaquette
double AverageSpatialPlaquette(vector<Eigen::Matrix3cd> const &links, vector<int> const &dims)
{
    // lattices sizes
    int Ns = dims[0];
    int Nt = dims[3];
    int NSite = cb(Ns) * Nt;

    // sum of spatial plaquettes
    double sum = 0.;
    // loop over spacetime coordinates and directions
    for (int x = 0; x < Ns; x++)
        for (int y = 0; y < Ns; y++)
            for (int z = 0; z < Ns; z++)
                for (int t = 0; t < Nt; t++)
                    for (int mu = 0; mu < 3; mu++)
                        for (int nu = mu + 1; nu < 3; nu++)
                        {
                            vector<int> coordinates = {x, y, z, t, mu};
                            sum += Plaquette(links, coordinates, dims, nu).trace().real();
                        }
    // returing average spatial plaquette
    return sum / 3. / (double)NSite;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// average temporal plaquette
double AverageTemporalPlaquette(vector<Eigen::Matrix3cd> const &links, vector<int> const &dims)
{
    // lattices sizes
    int Ns = dims[0];
    int Nt = dims[3];
    int NSite = cb(Ns) * Nt;

    // sum of temporal plaquettes
    double sum = 0.;
    // loop over spacetime coordinates and directions
    for (int x = 0; x < Ns; x++)
        for (int y = 0; y < Ns; y++)
            for (int z = 0; z < Ns; z++)
                for (int t = 0; t < Nt; t++)
                    for (int mu = 0; mu < 3; mu++)
                    {
                        int nu = 3;
                        vector<int> coordinates = {x, y, z, t, mu};
                        sum += Plaquette(links, coordinates, dims, nu).trace().real();
                    }
    // return average temporal plaquette
    return sum / 3. / (double)NSite;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// Polyakov loop
complex<double> PolyakovLoop(vector<Eigen::Matrix3cd> const &links, vector<int> const &dims, vector<int> const &spatialSite)
{
    // temporal size
    int Nt = dims[3];

    // product of links
    Eigen::Matrix3cd prod = Eigen::Matrix3cd::Identity();
    // loop over temporal coordinates at the given site
    for (int t = 0; t < Nt; t++)
    {
        // U_4(x, y, z, t)
        vector<int> coordinates = {spatialSite[0], spatialSite[1], spatialSite[2], t, 3};
        prod *= links[GetIndex(coordinates, dims)];
    }
    // return Polyakov loop
    return prod.trace();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// Polyakov loop withput taking the trace
Eigen::Matrix3cd PolyakovLoop_noTrace(vector<Eigen::Matrix3cd> const &links, vector<int> const &dims, vector<int> const &coordinates)
{
    // temporal size
    int Nt = dims[3];

    // product of links
    Eigen::Matrix3cd prod = ID3;
    // loop over temporal coordinates at the given site
    // starting loop at given link coordinates
    for (int t = coordinates[3]; t < Nt + coordinates[3]; t++)
    {
        // U_4(x, y, z, t)
        vector<int> linkCoordinates = {coordinates[0], coordinates[1], coordinates[2], t, 3};
        prod *= links[GetIndex(linkCoordinates, dims)];
    }
    // return Polyakov loop with taking the trace
    return prod;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// general Polyakov loop ~ loop in any dimension
complex<double> GeneralPolyakovLoop(vector<Eigen::Matrix3cd> const &links, vector<int> const &dims, vector<int> const &site, int const &dummyIndex)
{
    // product of links
    Eigen::Matrix3cd prod = Eigen::Matrix3cd::Identity();
    // loop over dummy coordinate at the given site
    int N = dims[dummyIndex];
    for (int dummy = 0; dummy < N; dummy++)
    {
        // U_dummy(...)
        vector<int> coordinates(5);
        coordinates[4] = dummyIndex;
        int siteCoord = 0;
        for (int i = 0; i < 4; i++)
        {
            if (i == dummyIndex)
                coordinates[i] = dummy;
            else
            {
                coordinates[i] = site[siteCoord];
                siteCoord++;
            }
        }
        prod *= links[GetIndex(coordinates, dims)];
    }
    // return Polyakov loop
    return prod.trace();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// average Polyakov loop
complex<double> AveragePolyakovLoop(vector<Eigen::Matrix3cd> const &links, vector<int> const &dims)
{
    // lattices sizes
    int Ns = dims[0];

    // sum of Polyakov loops
    complex<double> sum = 0.;
    // loop over spatial coordinates
    for (int x = 0; x < Ns; x++)
        for (int y = 0; y < Ns; y++)
            for (int z = 0; z < Ns; z++)
            {
                vector<int> spatialSite = {x, y, z};
                sum += PolyakovLoop(links, dims, spatialSite);
            }
    // return average Polyakov loop
    return sum / (double)cb(Ns);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// average spatial / space-like Polyakov loop
complex<double> AverageSpatialPolyakovLoop(vector<Eigen::Matrix3cd> const &links, vector<int> const &dims)
{
    // lattices sizes
    int Ns = dims[0];
    int Nt = dims[3];

    // sum of Polyakov loops
    complex<double> sum = 0.;
    // loop over spatial dimensions
    for (int d = 0; d < 3; d++)
    {
        // loop over non-dummy spatial and temporal coordinates
        for (int s1 = 0; s1 < Ns; s1++)
            for (int s2 = 0; s2 < Ns; s2++)
                for (int t = 0; t < Nt; t++)
                {
                    vector<int> site = {s1, s2, t};
                    sum += GeneralPolyakovLoop(links, dims, site, d);
                }
    }
    // return average spatial Polyakov loop
    return sum / (double)(sq(Ns) * Nt) / 3.;
}

#pragma warning(pop)