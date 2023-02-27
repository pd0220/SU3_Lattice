#pragma once

#pragma warning(push)
#pragma warning(disable : 4127)

// used headers and/or libraries
#include <random>

// custom
#include "auxiliary.hh"
#include "lattice_io.hh"
#include "pureGauge.hh"
#include "HDQCD.hh"
#include "observables.hh"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// _MAIN_ function
int main(int argc, char **argv)
{
    // if not enough arguments given
    if (argc < 6)
    {
        cout << "ERROR in main: not enough arguments given." << endl;
        exit(-1);
    }

    // lattice sizes
    const int Ns = atoi(argv[1]);
    const int Nt = atoi(argv[2]);
    // number of links for 4 spacetime dimensions
    const int NSite = cb(Ns) * Nt;
    const int NLink = NSite * 4;

    // coupling
    const double beta = atof(argv[3]);
    // coefficient of the rectangular plaquette
    const double c1 = atof(argv[4]);
    // coefficient of the squared plaquette
    const double c2 = atof(argv[5]);
    // coefficient of the plaquette from normalization: c0 + 8 c1 = 1
    const double c0 = 1. - 8. * c1;
    // spread parameter
    const double eps = atof(argv[6]);
    // number of MC sweeps
    const int T = atoi(argv[7]);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // random number generation
    std::random_device rd{};
    std::mt19937 gen(rd());
    // uniform random from [0, 1)
    std::uniform_real_distribution<double> uniformDistr(0., 1.);
    auto RandomNumber = [&]()
    {
        return uniformDistr(gen);
    };
    // uniform random from [0, Ns)
    std::uniform_int_distribution uniformNsDistr(0, Ns - 1);
    auto RandomNsNumber = [&]()
    {
        return uniformNsDistr(gen);
    };
    // uniform random from [0, Nt)
    std::uniform_int_distribution uniformNtDistr(0, Nt - 1);
    auto RandomNtNumber = [&]()
    {
        return uniformNtDistr(gen);
    };
    // unform random from [0, 4)
    std::uniform_int_distribution uniformDirDistr(0, 3);
    auto RandomDirNumber = [&]()
    {
        return uniformDirDistr(gen);
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // dimensions
    vector<int> dims = {Ns, Ns, Ns, Nt, 4};
    // initial coordinates
    vector<int> coordinates = {0, 0, 0, 0, 0};

    // cold start
    vector<Eigen::Matrix3cd> links(NLink, ID3);

    // configuration read test
    /*
    const char *p = "lat_12_4_181126_151350";
    vector<Eigen::Matrix3cd> links = read_lattice(const_cast<char *>(p));

    cout << AverageTemporalPlaquette(links, dims) << endl;
    cout << AverageSpatialPlaquette(links, dims) << endl;
    cout << AveragePolyakovLoop(links, dims) << endl;
    cout << AverageSpatialPolyakovLoop(links, dims) << endl;
    */

    // configuration save test
    /*
    const char *t = "save_lattice_TEST";
    save_lattice(const_cast<char *>(t), Ns, Nt, links);
    vector<Eigen::Matrix3cd> links2 = read_lattice(const_cast<char *>(t));

    cout << AverageTemporalPlaquette(links, dims) << endl;
    cout << AverageSpatialPlaquette(links, dims) << endl;
    cout << AveragePolyakovLoop(links, dims) << endl;
    cout << AverageSpatialPolyakovLoop(links, dims) << endl;
    */

    // stout smearing tests
    /*
    vector<Eigen::Matrix3cd> smearedLinks1 = StoutLattice(links, dims, 0.125);

    cout << AverageTemporalPlaquette(smearedLinks1, dims) << endl;
    cout << AverageSpatialPlaquette(smearedLinks1, dims) << endl;
    cout << AveragePolyakovLoop(smearedLinks1, dims) << endl;
    cout << AverageSpatialPolyakovLoop(smearedLinks1, dims) << endl;

    vector<Eigen::Matrix3cd> smearedLinks2 = StoutLattice(smearedLinks1, dims, 0.125);

    cout << AverageTemporalPlaquette(smearedLinks2, dims) << endl;
    cout << AverageSpatialPlaquette(smearedLinks2, dims) << endl;
    cout << AveragePolyakovLoop(smearedLinks2, dims) << endl;
    cout << AverageSpatialPolyakovLoop(smearedLinks2, dims) << endl;
    */

    // RUN simulation
    // SWEEPS
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
            Eigen::Matrix3cd X = UpdatingMatrix(eps, RandomNumber);

            // current link
            Eigen::Matrix3cd U = links[GetIndex(coordinates, dims)];

            // local change in the Wilson action
            double deltaGaugeAction = DeltaGaugeAction(beta, c0, c1, c2, links, coordinates, dims, U, X);
            // rate
            double rate = Rate(deltaGaugeAction);
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
                // links[m].row(2) = ((links[m].row(0).conjugate()).cross(links[m].row(1).conjugate())).conjugate();

            complex<double> avgPLoop = AveragePolyakovLoop(links, dims);
            complex<double> avgSPLoop = AverageSpatialPolyakovLoop(links, dims);
            cout << AveragePlaquette(links, dims) << " " << avgPLoop.real() << " " << avgPLoop.imag() << " " << avgSPLoop.real() << " " << avgSPLoop.imag() << endl;
        }
    }
}

#pragma warning(pop)