#pragma once
#pragma warning(push)
#pragma warning(disable : 4127)

// used headers and/or libraries
#include <iostream>
#include <fstream>

// custom
#include "auxiliary.hh"

using namespace std;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// reading links from standardized lattice file format
// type definition for the storage format ~ coefficients of Gell-Mann matrices
typedef struct
{
    float a[8];
} su3_matrix_io;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

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
    // 2 pi > 6 ~ always true
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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// convert SU(3) matric to Gell-Mann coefficients
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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// computes Tr(A lambda_n) where lambda_n is the nth Gell-Mann matrix
double dot_gellmann(Eigen::Matrix3cd const &A, int const &n)
{
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
        printf("No such Gell-Mann matrix exists: n= %d\b", n);
        exit(-1);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// convert Gell-Mann coefficients to SU(3) matrix
void convert_mat_ftoio(Eigen::Matrix3cd const &M, su3_matrix_io *result)
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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// reading lattice format (4 spacetime dimension)
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

    int NSite = x * y * z * t;
    int NLink = NSite * 4;

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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

// saving lattice format (4 spacetime dimensions)
void save_lattice(char *filename, int const &Ns, int const &Nt, vector<Eigen::Matrix3cd> const &links)
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

    int NSite = cb(Ns) * Nt;

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

#pragma warning(pop)