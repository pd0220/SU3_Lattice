#pragma once
#pragma warning(push)
#pragma warning(disable : 4127)

// used headers and/or libraries
// custom
#include "auxiliary.hh"

using namespace std;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

//
// NAVIGATON
// generalized row-major format is used for indexing the sites and links instead of a 5-dimensional array
//
// index function for easier understanding and moving around in the flattened array
// taking into account different possible lattice extents and different possible matrix sizes too
// (fixed to 4 spacetime dimensions)
//
// get index from coordinates
int GetIndex(vector<int> const &coordinates, vector<int> const &dims)
{
    // check if mu > 4
    if ((coordinates[4] > 3) || (coordinates[4] < 0))
    {
        cout << "ERROR in GetIndex: mu > 3 or mu < 0" << endl;
        exit(-1);
    }
    // using periodic boundary conditions
    vector<int> redefined_coordinates = coordinates;
    for (int i = 0; i < 4; i++)
    {
        // addig an extra term because C++ is stupid...
        redefined_coordinates[i] = (coordinates[i] + dims[i]) % dims[i];
    }
    // flattened index
    int flattenedIndex = redefined_coordinates[0];
    // strides ~ bigger steps
    int stride = 1;
    for (int i = 1; i < 5; i++)
    {
        stride *= dims[i - 1];
        flattenedIndex += redefined_coordinates[i] * stride;
    }

    // return flattened index
    return flattenedIndex;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// get coordinates and spacetime direction from index ~ (x, y, z, t, mu)
vector<int> GetCoordinates(int const &index, vector<int> const &dims)
{
    // array of coordinates and spacetime direction
    vector<int> coordinates = {0, 0, 0, 0, 0};
    // offset
    int offset = index;
    for (int i = 0; i < 5; i++)
    {
        coordinates[i] = offset % dims[i];
        offset /= dims[i];
    }

    // return coordinates and spacetime direction
    return coordinates;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

//
// MONTE CARLO
//

// upper staple (i.e. with given nu) at given site with given mu
// generate arguments of links
vector<Eigen::MatrixXi> StapleUpper_args(vector<int> const &coordinates, int const &nu)
{
    // given direction mu
    int mu = coordinates[4];

    // link coordinates
    // &
    // if 0 take adjoint
    // else do nothing
    Eigen::RowVectorXi coordinates_and_ifDagger(6);
    coordinates_and_ifDagger << coordinates[0], coordinates[1], coordinates[2], coordinates[3], coordinates[4], 1;

    // used link coordinate combinations in the rectangular staples
    Eigen::MatrixXi initArguments(coordinates_and_ifDagger.colwise().replicate(3));

    // used link coordinate combinations in the staple
    vector<Eigen::MatrixXi> arguments(1, initArguments);

    // arguments
    // U_nu(n + mu)
    arguments[0](0, mu) += 1;
    arguments[0](0, 4) = nu;
    // U_mu(n + nu)^dagger
    arguments[0](1, nu) += 1;
    arguments[0](1, 5) = 0;
    // U_nu(n)^dagger
    arguments[0](2, 4) = nu;
    arguments[0](2, 5) = 0;

    // return arguments
    return arguments;
    /*
    // construct staple
    // sum of rectangular staples
    Eigen::Matrix3cd staple = ID3;

    for (int link = 0; link < 3; link++)
    {
        // adjoint or not?
        if (arguments(link, 5) == 0)
            staple *= links[GetIndex(eigen_to_std(arguments.row(link)), dims)].adjoint();
        else
            staple *= links[GetIndex(eigen_to_std(arguments.row(link)), dims)];
    }
    // return sum of staples
    return staple;
    */
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// upper staple (i.e. with given nu) at given site with given mu
// generate arguments of links
vector<Eigen::MatrixXi> StapleLower_args(vector<int> const &coordinates, int const &nu)
{
    // given direction mu
    int mu = coordinates[4];

    // link coordinates
    // &
    // if 0 take adjoint
    // else do nothing
    Eigen::RowVectorXi coordinates_and_ifDagger(6);
    coordinates_and_ifDagger << coordinates[0], coordinates[1], coordinates[2], coordinates[3], coordinates[4], 1;

    // used link coordinate combinations in the rectangular staples
    Eigen::MatrixXi initArguments(coordinates_and_ifDagger.colwise().replicate(3));

    // used link coordinate combinations in the staple
    vector<Eigen::MatrixXi> arguments(1, initArguments);

    // arguments
    // U_nu(n + mu - nu)^dagger
    arguments[0](0, mu) += 1;
    arguments[0](0, nu) -= 1;
    arguments[0](0, 4) = nu;
    arguments[0](0, 5) = 0;
    // U_mu(n - nu)^dagger
    arguments[0](1, nu) -= 1;
    arguments[0](1, 5) = 0;
    // U_nu(n - nu)
    arguments[0](2, nu) -= 1;
    arguments[0](2, 4) = nu;

    // return arguments
    return arguments;
    /*
    // construct staple
    // sum of rectangular staples
    Eigen::Matrix3cd staple = ID3;

    for (int link = 0; link < 3; link++)
    {
        // adjoint or not?
        if (arguments(link, 5) == 0)
            staple *= links[GetIndex(eigen_to_std(arguments.row(link)), dims)].adjoint();
        else
            staple *= links[GetIndex(eigen_to_std(arguments.row(link)), dims)];
    }

    // return sum of staples
    return staple;
    */

    // return arguments
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// sum of upper and lower (i.e. with given nu) staples at given site with given mu ~ "A"
/*
Eigen::Matrix3cd StapleSum(vector<Eigen::Matrix3cd> const &links, vector<int> const &coordinates, vector<int> const &dims, int const &nu)
{
    // given mu direction
    int mu = coordinates[4];

    // used coordinate combinations in the staples

    // UPPER
    //
    // arguments
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

    // product
    Eigen::Matrix3cd First = links[GetIndex(n_plus_mu_at_nu, dims)] * links[GetIndex(n_plus_nu_at_mu, dims)].adjoint() * links[GetIndex(n_at_nu, dims)].adjoint();

    // LOWER
    //
    // arguments
    // U_nu(n + mu - nu)
    vector<int> n_plus_mu_minus_nu_at_nu = coordinates;
    n_plus_mu_minus_nu_at_nu[mu] += 1;
    n_plus_mu_minus_nu_at_nu[nu] -= 1;
    n_plus_mu_minus_nu_at_nu[4] = nu;
    // U_mu(n - nu)
    vector<int> n_minus_nu_at_mu = coordinates;
    n_minus_nu_at_mu[nu] -= 1;
    // U_nu(n - nu)
    vector<int> n_minus_nu_at_nu = coordinates;
    n_minus_nu_at_nu[nu] -= 1;
    n_minus_nu_at_nu[4] = nu;

    // product
    Eigen::Matrix3cd Second = links[GetIndex(n_plus_mu_minus_nu_at_nu, dims)].adjoint() * links[GetIndex(n_minus_nu_at_mu, dims)].adjoint() * links[GetIndex(n_minus_nu_at_nu, dims)];

    // return sum of staples
    return First + Second;
}
*/

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// sum of upper and lower (i.e. with given nu) RECTANGULAR staples at given site with given mu ~ "B"
// generate arguments of links
vector<Eigen::MatrixXi> StapleSumImproved_args(vector<int> const &coordinates, int const &nu)
{
    // given mu direction
    int mu = coordinates[4];

    // link coordinates
    // &
    // if 0 take adjoint
    // else do nothing
    Eigen::RowVectorXi coordinates_and_ifDagger(6);
    coordinates_and_ifDagger << coordinates[0], coordinates[1], coordinates[2], coordinates[3], coordinates[4], 1;

    // used link coordinate combinations in the rectangular staples
    Eigen::MatrixXi initArguments(coordinates_and_ifDagger.colwise().replicate(5));

    // argument table for the 5 different links in the 6 different staples
    vector<Eigen::MatrixXi> arguments(6, initArguments);

    // UPPER "standing"
    //
    // arguments
    // U_nu(n + mu)
    arguments[0](0, mu) += 1;
    arguments[0](0, 4) = nu;
    // U_nu(n + mu + nu)
    arguments[0](1, mu) += 1;
    arguments[0](1, nu) += 1;
    arguments[0](1, 4) = nu;
    // U_mu(n + 2 nu)^dagger
    arguments[0](2, nu) += 2;
    arguments[0](2, 5) = 0;
    // U_nu(n + nu)^dagger
    arguments[0](3, nu) += 1;
    arguments[0](3, 4) = nu;
    arguments[0](3, 5) = 0;
    // U_nu(n)^dagger
    arguments[0](4, 4) = nu;
    arguments[0](4, 5) = 0;

    // LOWER "standing"
    //
    // arguments
    // U_nu(n + mu - nu)^dagger
    arguments[1](0, mu) += 1;
    arguments[1](0, nu) -= 1;
    arguments[1](0, 4) = nu;
    arguments[1](0, 5) = 0;
    // U_nu(n + mu - 2 nu)^dagger
    arguments[1](1, mu) += 1;
    arguments[1](1, nu) -= 2;
    arguments[1](1, 4) = nu;
    arguments[1](1, 5) = 0;
    // U_mu(n - 2 nu)^dagger
    arguments[1](2, nu) -= 2;
    arguments[1](2, 5) = 0;
    // U_nu(n - 2 nu)
    arguments[1](3, nu) -= 2;
    arguments[1](3, 4) = nu;
    // U_nu(n - nu)
    arguments[1](4, nu) -= 1;
    arguments[1](4, 4) = nu;

    // UPPER "lying" to RIGHT
    //
    // arguments
    // U_mu(n + mu)
    arguments[2](0, mu) += 1;
    // U_nu(n + 2 mu)
    arguments[2](1, mu) += 2;
    arguments[2](1, 4) = nu;
    // U_mu(n + mu + nu)^dagger
    arguments[2](2, mu) += 1;
    arguments[2](2, nu) += 1;
    arguments[2](2, 5) = 0;
    // U_mu(n + nu)^dagger
    arguments[2](3, nu) += 1;
    arguments[2](3, 5) = 0;
    // U_nu(n)^dagger
    arguments[2](4, 4) = nu;
    arguments[2](4, 5) = 0;

    // UPPER "lying" to LEFT
    //
    // arguments
    // U_nu(n + mu)
    arguments[3](0, mu) += 1;
    arguments[3](0, 4) = nu;
    // U_mu(n + nu)^dagger
    arguments[3](1, nu) += 1;
    arguments[3](1, 5) = 0;
    // U_mu(n - mu + nu)^dagger
    arguments[3](2, mu) -= 1;
    arguments[3](2, nu) += 1;
    arguments[3](2, 5) = 0;
    // U_nu(n - mu)^dagger
    arguments[3](3, mu) -= 1;
    arguments[3](3, 4) = nu;
    arguments[3](3, 5) = 0;
    // U_mu(n - mu)
    arguments[3](4, mu) -= 1;

    // LOWER "lying" to RIGHT
    //
    // arguments
    // U_mu(n + mu)
    arguments[4](0, mu) += 1;
    // U_nu(n + 2 mu - nu)^dagger
    arguments[4](1, mu) += 2;
    arguments[4](1, nu) -= 1;
    arguments[4](1, 4) = nu;
    arguments[4](1, 5) = 0;
    // U_mu(n + mu - nu)^dagger
    arguments[4](2, mu) += 1;
    arguments[4](2, nu) -= 1;
    arguments[4](2, 5) = 0;
    // U_mu(n - nu)^dagger
    arguments[4](3, nu) -= 1;
    arguments[4](3, 5) = 0;
    // U_nu(n - nu)
    arguments[4](4, nu) -= 1;
    arguments[4](4, 4) = nu;

    // LOWER "lying" to LEFT
    //
    // arguments
    // U_nu(n + mu - nu)^dagger
    arguments[5](0, mu) += 1;
    arguments[5](0, nu) -= 1;
    arguments[5](0, 4) = nu;
    arguments[5](0, 5) = 0;
    // U_mu(n - nu)^dagger
    arguments[5](1, nu) -= 1;
    arguments[5](1, 5) = 0;
    // U_mu(n - mu - nu)^dagger
    arguments[5](2, mu) -= 1;
    arguments[5](2, nu) -= 1;
    arguments[5](2, 5) = 0;
    // U_nu(n - mu - nu)
    arguments[5](3, mu) -= 1;
    arguments[5](3, nu) -= 1;
    arguments[5](3, 4) = nu;
    // U_mu(n - mu)
    arguments[5](4, mu) -= 1;

    // return arguments;
    return arguments;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// sum of of RECTANGULAR staples
/*
Eigen::Matrix3cd StapleSumImproved(vector<Eigen::Matrix3cd> const &links, vector<int> const &coordinates, vector<int> const &dims, int const &nu)
{
    // generate arguments
    vector<Eigen::MatrixXi> arguments = StapleSumImproved_args(coordinates, nu);
    // sum of rectangular staples
    Eigen::Matrix3cd stapleSum = Eigen::Matrix3cd::Zero();
    Eigen::Matrix3cd tmpMat = ID3;

    for (int staple = 0; staple < 6; staple++)
    {
        // construct staple
        for (int link = 0; link < 5; link++)
        {
            // adjoint or not?
            if (arguments[staple](link, 5) == 0)
                tmpMat *= links[GetIndex(eigen_to_std(arguments[staple].row(link)), dims)].adjoint();
            else
                tmpMat *= links[GetIndex(eigen_to_std(arguments[staple].row(link)), dims)];
        }
        // add up staples
        stapleSum += tmpMat;
        tmpMat = ID3;
    }
    // return sum of staples
    return stapleSum;
}
*/

// sum of general staples from sets of link arguments
Eigen::Matrix3cd StapleSum(vector<Eigen::Matrix3cd> const &links, vector<Eigen::MatrixXi> const &arguments, vector<int> const &dims)
{
    // initialize sum of staples
    Eigen::Matrix3cd stapleSum = Eigen::Matrix3cd::Zero();
    // auxiliary matrix
    Eigen::Matrix3cd tmpMat = ID3;
    // number of staples in the sum
    int NumStaples = static_cast<int>(arguments.size());
    // number of links in a staple
    int NumLinks = arguments[0].rows();
    // loop over staples
    for (int staple = 0; staple < NumStaples; staple++)
    {
        // construct staple
        for (int link = 0; link < NumLinks; link++)
        {
            // adjoint or not?
            if (arguments[staple](link, 5) == 0)
                tmpMat *= links[GetIndex(eigen_to_std(arguments[staple].row(link)), dims)].adjoint();
            else
                tmpMat *= links[GetIndex(eigen_to_std(arguments[staple].row(link)), dims)];
        }
        // add up staples
        stapleSum += tmpMat;
        // reset auxiliary matrix
        tmpMat = ID3;
    }
    // return sum of staples
    return stapleSum;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// generating SU(2) matrix
Eigen::Matrix2cd SU2(double const &eps, function<double()> const &randFunc)
{
    // generating four random numbers from [-1/2, 1/2]
    vector<double> r(4, 0);
    for (int i = 0; i < 4; i++)
        r[i] = randFunc() - 0.5;

    // length of 3-vector
    double norm = sqrt(sq(r[0]) + sq(r[1]) + sq(r[2]));

    // generate coefficients
    vector<double> x(4, 0);
    for (int i = 0; i < 3; i++)
        x[i] = eps * r[i] / norm;
    x[3] = ((r[3] > 0.) - (r[3] < 0)) * sqrt(1. - sq(eps));

    // parameterizing SU(2) matrices
    Eigen::Matrix2cd X2 = x[3] * ID2;
    for (int i = 0; i < 3; i++)
        X2 += I * Pauli[i] * x[i];

    // return the SU(2) matrix
    return X2;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// candidate link generation through updating matrix X
Eigen::Matrix3cd UpdatingMatrix(double const &eps, function<double()> randFunc)
{
    // generating three SU(2) matrices
    vector<Eigen::Matrix2cd> SU2Matrices(3, Eigen::Matrix2cd::Zero());
    for (int i = 0; i < 3; i++)
        SU2Matrices[i] = SU2(eps, randFunc);

    // generating three 3x3 matrices ~ SU(2) subgropups method
    Eigen::Matrix3cd R{{SU2Matrices[0](0, 0), SU2Matrices[0](0, 1), complex(0., 0.)},
                       {SU2Matrices[0](1, 0), SU2Matrices[0](1, 1), complex(0., 0.)},
                       {complex(0., 0.), complex(0., 0.), complex(1., 0.)}};

    Eigen::Matrix3cd S{{SU2Matrices[1](0, 0), complex(0., 0.), SU2Matrices[1](0, 1)},
                       {complex(0., 0.), complex(1., 0.), complex(0., 0.)},
                       {SU2Matrices[1](1, 0), complex(0., 0.), SU2Matrices[1](1, 1)}};

    Eigen::Matrix3cd T{{complex(1., 0.), complex(0., 0.), complex(0., 0.)},
                       {complex(0., 0.), SU2Matrices[2](0, 0), SU2Matrices[2](0, 1)},
                       {complex(0., 0.), SU2Matrices[2](1, 0), SU2Matrices[2](1, 1)}};

    // return the X matrix
    // TODO
    // choose randomly between X and its adjoint
    if (randFunc() > 0.5)
        return R * S * T;
    else
        return (R * S * T).adjoint();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// local change in Wilson action (single site)
// possible TO-DO
// rewrite to use U' = X * U instead of U and X separately
double DeltaGaugeAction(double const &beta, double const &c0, double const &c1, double const &c2, vector<Eigen::Matrix3cd> const &links, vector<int> const &coordinates, vector<int> const &dims, Eigen::Matrix3cd const &U, Eigen::Matrix3cd const &X)
{
    // Lorentz index of the given link
    int mu = coordinates[4];

    // sum of staples
    // Eigen::Matrix3cd SumOfStaples = Eigen::Matrix3cd::Zero();
    Eigen::Matrix3cd A = Eigen::Matrix3cd::Zero();
    Eigen::Matrix3cd B = Eigen::Matrix3cd::Zero();
    Eigen::Matrix3cd AUpper = Eigen::Matrix3cd::Zero();
    Eigen::Matrix3cd ALower = Eigen::Matrix3cd::Zero();

    for (int nu = 0; nu < 4; nu++)
    {
        if (nu == mu)
            continue;
        else
        {
            // A += StapleSum(links, coordinates, dims, nu);
            AUpper += StapleSum(links, StapleUpper_args(coordinates, nu), dims);
            ALower += StapleSum(links, StapleLower_args(coordinates, nu), dims);
            B += StapleSum(links, StapleSumImproved_args(coordinates, nu), dims);
        }
    }
    // sum of lower and upper plaquette staples
    A = AUpper + ALower;

    // return local change in Wilson action
    // return -beta / 3. * ((X - ID3) * U * A).real().trace();
    // return local change in the improved Wilson action
    // return -beta / 9. * ((X - ID3) * U * (5 * A - 0.25 * B)).real().trace();
    // return -beta / 3. * ((X - ID3) * U * (c0 * A + c1 * B)).real().trace();

    // differences upon proposal in the plaquette and the improved rectangular terms
    Eigen::Matrix3cd diffLin = (X - ID3) * U * (c0 * A + c1 * B);

    // differences upon proposal in the squared plaquette term
    // M_i = U * A_i
    Eigen::Matrix3cd MUpper = U * AUpper;
    Eigen::Matrix3cd MLower = U * ALower;
    // computing the difference
    Eigen::Matrix3cd diffSq = c2 * (((X * MUpper * X) - MUpper) * MUpper + ((X * MLower * X) - MLower) * MLower);

    // return local change in action
    return -beta / 3. * (diffLin + diffSq).real().trace();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// acceptance rate for the Metropolis algorithm
double Rate(double const &deltaAction)
{
    if (deltaAction <= 0)
        return 1.;
    else
        return std::exp(-deltaAction);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

//
// STOUT SMEARING
// Phys.Rev. D69 (2004) 054501
//

// arguments for the staple sum involved in the smearing step
vector<Eigen::MatrixXi> StapleStout_args(vector<int> const &coordinates, int const &nu)
{
    // given direction mu
    int mu = coordinates[4];

    // link coordinates
    // &
    // if 0 take adjoint
    // else do nothing
    Eigen::RowVectorXi coordinates_and_ifDagger(6);
    coordinates_and_ifDagger << coordinates[0], coordinates[1], coordinates[2], coordinates[3], coordinates[4], 1;

    // used link coordinate combinations in the rectangular staples
    Eigen::MatrixXi initArguments(coordinates_and_ifDagger.colwise().replicate(3));

    // used link coordinate combinations in the staple
    vector<Eigen::MatrixXi> arguments(2, initArguments);

    // UPPER
    // arguments
    // U_nu(n)
    arguments[0](0, 4) = nu;
    // U_mu(n + nu)
    arguments[0](1, nu) += 1;
    // U_nu(n + mu)^dagger
    arguments[0](2, mu) += 1;
    arguments[0](2, 4) = nu;
    arguments[0](2, 5) = 0;

    // LOWER
    // arguments
    // U_nu(n - nu)^dagger
    arguments[1](0, nu) -= 1;
    arguments[1](0, 4) = nu;
    arguments[1](0, 5) = 0;
    // U_mu(n - nu)
    arguments[1](1, nu) -= 1;
    // U_nu(n - nu + mu)
    arguments[1](2, mu) += 1;
    arguments[1](2, nu) -= 1;
    arguments[1](2, 4) = nu;

    // return arguments
    return arguments;
}

// performing a single smearing step on a single link
// TO-DO ~ handling exceptions...
// isotropic: rho_{mu, nu} = rho
Eigen::Matrix3cd StoutLink(vector<Eigen::Matrix3cd> const &links, vector<int> const &coordinates, vector<int> const &dims, double const &rho)
{
    // Lorentz index of the given link
    int mu = coordinates[4];

    // construct staple sum ~ C_{mu, nu}
    Eigen::Matrix3cd C = Eigen::Matrix3cd::Zero();
    for (int nu = 0; nu < 4; nu++)
    {
        if (nu == mu)
            continue;
        else
            C += StapleSum(links, StapleStout_args(coordinates, nu), dims);
    }
    C *= rho;

    // link in question
    Eigen::Matrix3cd U = links[GetIndex(coordinates, dims)];
    // construct Omega = C * U^dagger
    Eigen::Matrix3cd Omega = C * U.adjoint();
    // Omega^dagger - Omega
    Eigen::Matrix3cd OmegaDiff = Omega.adjoint() - Omega;

    // construct Q
    Eigen::Matrix3cd Q = I / 2. * (OmegaDiff - ID3 * OmegaDiff.trace() / 3.);

    // return smeared link
    return (I * Q).exp() * U;
}

// stout smearing the whole lattice
vector<Eigen::Matrix3cd> StoutLattice(vector<Eigen::Matrix3cd> const &links, vector<int> const &dims, double const &rho)
{
    // number of links
    int NumLinks = static_cast<int>(links.size());
    // stout smeared lattice
    vector<Eigen::Matrix3cd> stoutLinks(links.size(), ID3);

    // loop over original links
    for (int link = 0; link < NumLinks; link++)
    {
        stoutLinks[link] = StoutLink(links, GetCoordinates(link, dims), dims, rho);
    }

    // return smeared lattice
    return stoutLinks;
}

#pragma warning(pop)