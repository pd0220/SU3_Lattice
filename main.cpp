#pragma warning(push)
#pragma warning(disable : 4127)

// used headers and/or libraries
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>
#include <complex>

// Eigen
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// imaginary unit
const complex<double> I(0., 1.);
// square root of 3
const double ONE_OVER_SQRT3 = 1. / sqrt(3.);
// complex 2x2 and 3x3 identity matrices
// pi
const double PI = 3.14159265359;
const Eigen::Matrix2cd ID2 = Eigen::Matrix2cd::Identity();
const Eigen::Matrix3cd ID3 = Eigen::Matrix3cd::Identity();

// Pauli matrices
const Eigen::Matrix2cd sigma1{{0, 1}, {1, 0}};
const Eigen::Matrix2cd sigma2{{0, -I}, {I, 0}};
const Eigen::Matrix2cd sigma3{{1, 0}, {0, -1}};
const vector<Eigen::Matrix2cd> Pauli = {sigma1, sigma2, sigma3};

// auxiliary lambda(s)
// general square function
auto sq = [](auto const &x)
{
    return x * x;
};

// general cubic function
auto cb = [](auto const &x)
{
    return x * x * x;
};

// general quartic function
auto quart = [](auto const &x)
{
    return sq(sq(x));
};

// lattice size
const int Ns = 6;
const int Nt = 4;
// number of links(?) for 4 spacetime dimensions
const int NSite = cb(Ns) * Nt;
const int NLink = NSite * 4;

// Eigen vector to std vector
vector<int> eigen_to_std(Eigen::VectorXi const &Eigen_v)
{
    vector<int> std_v(Eigen_v.data(), Eigen_v.data() + Eigen_v.size());
    return std_v;
}

// random number generation
std::random_device rd{};
std::mt19937 gen(rd());
// uniform random from [0, 1)
std::uniform_real_distribution<double> uniformDistr(0., 1.);
double RandomNumber()
{
    return uniformDistr(gen);
}
// uniform random from [0, Ns)
std::uniform_int_distribution uniformNsDistr(0, Ns - 1);
int RandomNsNumber()
{
    return uniformNsDistr(gen);
}
// uniform random from [0, Nt)
std::uniform_int_distribution uniformNtDistr(0, Nt - 1);
int RandomNtNumber()
{
    return uniformNtDistr(gen);
}
// unform random from [0, 4)
std::uniform_int_distribution uniformDirDistr(0, 3);
int RandomDirNumber()
{
    return uniformDirDistr(gen);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

    // retur coordinates and spacetime direction
    return coordinates;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// trace free logarithm of a (unitary) matrix
// explicitly using that we have SU(3) matrices
// [arxiv: 0709.4110] (A.7)
Eigen::Matrix3cd TraceFreeLogM(Eigen::Matrix3cd const &U, double const &ImTrace)
{
    // solving the eigenvalue problem
    Eigen::ComplexEigenSolver<Eigen::Matrix3cd> CES;
    CES.compute(U);
    Eigen::Vector3cd eVals = CES.eigenvalues();
    Eigen::Matrix3cd eVecs = CES.eigenvectors();
    // taking the logarithm of each eigenvalue
    Eigen::Vector3cd log_eVals = Eigen::log(eVals.array());
    // performing redefinition of the largest eigenvalue through the subtraction
    double t = ImTrace;

    // initialize logarithm
    Eigen::Matrix3cd logM = Eigen::Matrix3cd::Zero();
    // 2pi > 6 always
    while (abs(t) > 6.)
    {
        // find the eigenvalue with the largest imaginary part
        Eigen::Index maxRow, maxCol;
        log_eVals.imag().maxCoeff(&maxRow, &maxCol);
        // acting according to sign of ImTrlogM
        if (t > 0)
            log_eVals(maxRow) -= 2. * PI * I;
        else
            log_eVals(maxRow) += 2. * PI * I;

        // computing matrix logarithm from projectors
        logM = Eigen::Matrix3cd::Zero();
        for (int i = 0; i < 3; i++)
        {
            logM += log_eVals(i) * eVecs.col(i) * eVecs.col(i).adjoint();
        }

        // update
        t = logM.trace().imag();
    }

    // return the trace free logarithm of the given matrix
    return logM;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// LEGACY
// reading links from standardized lattice file format
// type definition for the storage format ~ coefficients of Gell-Mann matrices
typedef struct
{
    float a[8];
} su3_matrix_io;

Eigen::Matrix3cd convert_mat_iotof(su3_matrix_io *alpha)
{
    Eigen::Matrix3cd A = Eigen::Matrix3cd::Zero();

    A(0, 0) = 0. + I * ((double)alpha->a[7] * ONE_OVER_SQRT3 + (double)alpha->a[2]) / 2.;
    A(0, 1) = (double)alpha->a[1] / 2. + I * (double)alpha->a[0] / 2.;
    A(0, 2) = (double)alpha->a[4] / 2. + I * (double)alpha->a[3] / 2.;
    A(1, 0) = -(double)alpha->a[1] / 2. + I * (double)alpha->a[0] / 2.;
    A(1, 1) = 0. + I * (-(double)alpha->a[2] + (double)alpha->a[7] * ONE_OVER_SQRT3) / 2.;
    A(1, 2) = (double)alpha->a[6] / 2. + I * (double)alpha->a[5] / 2.;
    A(2, 0) = -(double)alpha->a[4] / 2. + I * (double)alpha->a[3] / 2.;
    A(2, 1) = -(double)alpha->a[6] / 2. + I * (double)alpha->a[5] / 2.;
    A(2, 2) = 0. + I * (-(double)alpha->a[7] * ONE_OVER_SQRT3);

    return A.exp();
}

double dot_gellmann(Eigen::Matrix3cd const &A, int n)
{
    // Computes Tr(A lambda_n) where lambda_n is the nth Gell-Mann matrix
    switch (n)
    {
    case 1:
        return (A(0, 1) + A(1, 0)).real();
        break;
    case 2:
        return (I * A(0, 1) - I * A(1, 0)).real();
        break;
    case 3:
        return (A(0, 0) - A(1, 1)).real();
        break;
    case 4:
        return (A(0, 2) + A(2, 0)).real();
        break;
    case 5:
        return (I * A(0, 2) - I * A(2, 0)).real();
        break;
    case 6:
        return (A(1, 2) + A(2, 1)).real();
        break;
    case 7:
        return (I * A(1, 2) - I * A(2, 1)).real();
        break;
    case 8:
        return (ONE_OVER_SQRT3 * (A(0, 0) + A(1, 1) - 2. * A(2, 2))).real();
        break;
    default:
        printf("No such Gell-Mann matrix exists; n= %d\b", n);
        exit(-1);
    }
}

void convert_mat_ftoio(Eigen::Matrix3cd M, su3_matrix_io *result)
{
    // built in logarithm
    Eigen::Matrix3cd logM = M.log();
    // if ImTrlogM = 2pi --> "shift Riemann sheet"
    if (abs(logM.trace().imag()) > 6.)
        logM = TraceFreeLogM(M, logM.trace().imag());
    Eigen::Matrix3cd L = -I * logM;
    for (int i = 0; i < 8; i++)
    {
        result->a[i] = static_cast<float>(dot_gellmann(L, i + 1));
    }
}

vector<Eigen::Matrix3cd> read_lattice(char *filename)
{
    FILE *f = NULL;

    if ((f = fopen(filename, "rb")) == 0)
    {
        printf("ERROR read_lattice: Unable to open file %s\n", filename);
        fflush(0);
        exit(1);
    }

    int x, y, z, t;

    fread(&x, sizeof(int), 1, f);
    fread(&y, sizeof(int), 1, f);
    fread(&z, sizeof(int), 1, f);
    fread(&t, sizeof(int), 1, f);

    printf("Lattice size= %dX%dX%dX%d\n", x, y, z, t);

    vector<Eigen::Matrix3cd> links(NLink, ID3);

    su3_matrix_io *complink = (su3_matrix_io *)malloc(NLink * sizeof(su3_matrix_io));

    fread(complink, sizeof(su3_matrix_io), NLink, f);

    for (int isite = 0; isite < NSite; isite++)
    {
        for (int idir = 0; idir < 4; idir++)
        {
            links[NSite * idir + isite] = convert_mat_iotof(complink + 4 * isite + idir);
        }
    }
    fclose(f);
    free(complink);
    return links;
}

void save_lattice(char *filename, vector<Eigen::Matrix3cd> links)
{
    printf("save_lattice: %s\n", filename);
    FILE *f = NULL;

    if ((f = fopen(filename, "wb")) == 0)
    {
        printf("ERROR save_lattice: Unable to open file %s\n", filename);
        fflush(0);
        exit(1);
    }

    fwrite(&Ns, sizeof(int), 1, f);
    fwrite(&Ns, sizeof(int), 1, f);
    fwrite(&Ns, sizeof(int), 1, f);
    fwrite(&Nt, sizeof(int), 1, f);

    su3_matrix_io complink;

    for (int isite = 0; isite < NSite; isite++)
    {
        for (int idir = 0; idir < 4; idir++)
        {
            convert_mat_ftoio(links[NSite * idir + isite], &complink);
            fwrite(&complink, sizeof(su3_matrix_io), 1, f);

            // free(complink);
        }
    }
    fclose(f);
    // free(complink);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// sum of upper and lower (i.e. with given nu) staples at given site with given mu ~ "A"
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
    // U(n - nu)
    vector<int> n_minus_nu_at_nu = coordinates;
    n_minus_nu_at_nu[nu] -= 1;
    n_minus_nu_at_nu[4] = nu;

    // product
    Eigen::Matrix3cd Second = links[GetIndex(n_plus_mu_minus_nu_at_nu, dims)].adjoint() * links[GetIndex(n_minus_nu_at_mu, dims)].adjoint() * links[GetIndex(n_minus_nu_at_nu, dims)];

    // return sum of staples
    return First + Second;
}

// sum of upper and lower (i.e. with given nu) RECTANGULAR staples at given site with given mu ~ "B"
// generate arguments of links
vector<Eigen::MatrixXi> StapleSumImproved_ARGUMENTS(vector<int> const &coordinates, int const &nu)
{
    // given mu direction
    int mu = coordinates[4];

    // link coordinates
    // &
    // if 0 take adjoint
    // else do not
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
    // U_mu(n - mu)^dagger
    arguments[5](1, mu) -= 1;
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

// sum of of RECTANGULAR staples
Eigen::Matrix3cd StapleSumImproved(vector<Eigen::Matrix3cd> const &links, vector<int> const &coordinates, vector<int> const &dims, int const &nu)
{
    // generate arguments
    vector<Eigen::MatrixXi> arguments = StapleSumImproved_ARGUMENTS(coordinates, nu);
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// generating SU(2) matrix
Eigen::Matrix2cd SU2(double const &eps)
{
    // generating four random numbers from [-1/2, 1/2]
    vector<double> r(4, 0);
    for (int i = 0; i < 4; i++)
        r[i] = RandomNumber() - 0.5;

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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// candidate link generation through updating matrix X
Eigen::Matrix3cd UpdatingMatrix(double const &eps)
{
    // generating three SU(2) matrices
    vector<Eigen::Matrix2cd> SU2Matrices(3, Eigen::Matrix2cd::Zero());
    for (int i = 0; i < 3; i++)
        SU2Matrices[i] = SU2(eps);

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
    if (RandomNumber() > 0.5)
        return R * S * T;
    else
        return (R * S * T).adjoint();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// local change in Wilson action (single site?)
double DeltaAction(double const &beta, double const &c0, double const &c1, vector<Eigen::Matrix3cd> const &links, vector<int> const &coordinates, vector<int> const &dims, Eigen::Matrix3cd const &U, Eigen::Matrix3cd const &X)
{
    // Lorentz index of the given link
    int mu = coordinates[4];
    // sum of staples
    // Eigen::Matrix3cd SumOfStaples = Eigen::Matrix3cd::Zero();
    Eigen::Matrix3cd A = Eigen::Matrix3cd::Zero();
    Eigen::Matrix3cd B = Eigen::Matrix3cd::Zero();

    for (int nu = 0; nu < 4; nu++)
    {
        if (nu == mu)
            continue;
        else
        {
            A += StapleSum(links, coordinates, dims, nu);
            B += StapleSumImproved(links, coordinates, dims, nu);
        }
    }

    // return local change Wilson action
    // return -beta / 3. * ((X - ID3) * U * A).real().trace();
    // return local change in the improved Wilson action
    // return -beta / 9. * ((X - ID3) * U * (5 * A - 0.25 * B)).real().trace();
    return -beta / 3. * ((X - ID3) * U * (c0 * A + c1 * B)).real().trace();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

// average plaquette
double AveragePlaquette(vector<Eigen::Matrix3cd> const &links, vector<int> const &dims)
{
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

// average spatial plaquette
double AverageSpatialPlaquette(vector<Eigen::Matrix3cd> const &links, vector<int> const &dims)
{
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

// average temporal plaquette
double AverageTemporalPlaquette(vector<Eigen::Matrix3cd> const &links, vector<int> const &dims)
{
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

// Polyakov loop
complex<double> PolyakovLoop(vector<Eigen::Matrix3cd> const &links, vector<int> const &dims, vector<int> const &spatialSite)
{
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

// average Polyakov loop
complex<double> AveragePolyakovLoop(vector<Eigen::Matrix3cd> const &links, vector<int> const &dims)
{
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

// average spatial / space-like Polyakov loop
complex<double> AverageSpatialPolyakovLoop(vector<Eigen::Matrix3cd> const &links, vector<int> const &dims)
{
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// acceptance rate for the Metropolis algorithm
double Rate(double const &deltaAction)
{
    if (deltaAction <= 0)
        return 1.;
    else
        return std::exp(-deltaAction);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// _MAIN_ function
int main(int argc, char **argv)
{
    // if not enough arguments given
    if (argc < 4)
    {
        cout << "ERROR in main: not enough arguments given." << endl;
        exit(-1);
    }
    // coupling
    double beta = atof(argv[1]);
    // coefficient of the rectangular plaquette
    double c1 = atof(argv[2]);
    // coefficient of the plaquette from normalization: c0 + 8 c1 = 1
    double c0 = 1. - 8. * c1;
    // spread parameter
    double eps = atof(argv[3]);
    // dimensions
    vector<int> dims = {Ns, Ns, Ns, Nt, 4};
    // initial coordinates
    vector<int> coordinates = {0, 0, 0, 0, 0};

    // cold start
    vector<Eigen::Matrix3cd> links(NLink, ID3);

    // testing
    // const char *p = "lat_12_4_181126_151350";
    // vector<Eigen::Matrix3cd> links = read_lattice(const_cast<char *>(p));

    // cout << AverageTemporalPlaquette(links, dims) << endl;
    // cout << AverageSpatialPlaquette(links, dims) << endl;
    // cout << AveragePolyakovLoop(links, dims) << endl;
    // cout << AverageSpatialPolyakovLoop(links, dims) << endl;

    // const char *t = "save_lattice_TEST";
    // save_lattice(const_cast<char *>(t), links);
    // vector<Eigen::Matrix3cd> links2 = read_lattice(const_cast<char *>(t));

    // RUN simulation
    // SWEEPS
    int T = (int)10000;
    int tau = 5;
    for (int t = 0; t < T; t++)
    {
        // single sweep via Metropolis steps
        for (int step = 0; step < NLink; step++)
        {
            // Metropolis step
            //
            // choosing a site and a direction randomly
            for (int coord = 0; coord < 3; coord++)
                coordinates[coord] = RandomNsNumber();
            coordinates[3] = RandomNtNumber();
            coordinates[4] = RandomDirNumber();

            // updating X
            Eigen::Matrix3cd X = UpdatingMatrix(eps);

            // current link
            Eigen::Matrix3cd U = links[GetIndex(coordinates, dims)];

            // local change in the Wilson action
            double deltaAction = DeltaAction(beta, c0, c1, links, coordinates, dims, U, X);
            // rate
            double rate = Rate(deltaAction);
            // random number from [0, 1)
            double r = RandomNumber();

            // decide if new link is accepted or not
            if (r < rate)
                links[GetIndex(coordinates, dims)] = X * U;
        }

        // MEASUREMENTS
        if ((t % tau) == 0)
        {
            // PROJECTING BACK TO SU(3)
            // for (int m = 0; m < NLink; m++)
            //    links[m].row(2) = (links[m].row(0).conjugate()).cross(links[m].row(1).conjugate());
            // measuring the average plaquette
            complex<double> avgPLoop = AveragePolyakovLoop(links, dims);
            complex<double> avgSPLoop = AverageSpatialPolyakovLoop(links, dims);
            cout << AveragePlaquette(links, dims) << " " << avgPLoop.real() << " " << avgPLoop.imag() << " " << avgSPLoop.real() << " " << avgSPLoop.imag() << endl;
        }
    }
}

#pragma warning(pop)