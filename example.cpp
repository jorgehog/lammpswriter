#include <iostream>
#include <sys/types.h>
#include <armadillo>
#include <mpi.h>

//Defining NDEBUG removes all checks (see list in the main function).
//#define NDEBUG

//Defining LAMMPSWRITER_USE_MPI enabled MPI support (_n, where n = rank, appended to file name)
//#define LAMMPSWRITER_USE_MPI

//Defining LAMMSWRITER_LOW_MEMORY enabled such that multiple data sets never gets stored simultaneously on a single node.
//This yields a slower file dump, but ensures stable and pain free file dump regardless of memory available and number of nodes active.
//#define LAMMPSWRITER_LOW_MEMORY

#include "lammpswriter/lammpswriter.h"

using namespace arma;
using namespace std;


void dumpLammps(const uint frameNumber,
                const uint nParticles,
                const uint nParticleProperties,
                const string path,
                const string name,
                const mat &pos,
                const mat &vel,
                const vec &type);

int main()
{
    /*
     * If NDEBUG is not defined, checks will be made to ensure that:
     *
     * 1. The set system size is consistent.
     * 2. The number of saved properties matches that specified.
     * 3. The file is properly handled between each frame.
     *
     */

    int rank = 0;
#ifdef LAMMPSWRITER_USE_MPI
    int nProcs;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srand(21321321 + rank); //Set a unique seed per process

    //Set the MPI rank and the master rank
    lammpswriter::setMPIRank(rank, nProcs);
#endif

    //Set the system size once and for all.
    //Can also set a start position and shear.
    lammpswriter::setSystemSize(1, 1, 1);

    uint nParticles = 10;
    nParticles += rank;

    //We store positions, velocities, and a type
    mat pos(nParticles, 3);
    mat vel(nParticles, 3);
    vec type = {1.0, 2.0};


    uint nParticleProperties = 7;

    string path = "/tmp";
    string name = "particleData";


    //Create an object for writing (used only for alternative one)
    lammpswriter writer(nParticleProperties, path, name);


    uint nCycles = 10;

    wall_clock timer;
    timer.tic();

    //First alternative usage
    for (uint c = 0; c < nCycles; ++c)
    {

        pos.randu();
        vel.randn();

        //for sequential usage, we initialize a new file specifying.
        writer.initializeNewFile(c, nParticles);

        for (uint i = 0; i < nParticles; ++i)
        {
            //Store the data we specified.
            writer << type[i%2]
                   << pos(i, 0)
                   << pos(i, 1)
                   << pos(i, 2)
                   << vel(i, 0)
                   << vel(i, 1)
                   << vel(i, 2);
        }

        //closes the file.
        writer.finalize();
    }

    if (rank == 0)
    {
        cout << "dump took " << timer.toc() << "seconds." << endl;
        timer.tic();
    }

    //Second alternative usage
    for (uint c = nCycles; c < 2*nCycles; ++c)
    {
        pos.randu();
        vel.randn();

        dumpLammps(c,
                   nParticles,
                   nParticleProperties,
                   path,
                   name,
                   pos,
                   vel,
                   type);
    }

    if (rank == 0)
    {
        cout << "dump took " << timer.toc() << "seconds." << endl;
        timer.tic();
    }

#ifdef LAMMPSWRITER_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}


void dumpLammps(const uint frameNumber,
                const uint nParticles,
                const uint nParticleProperties,
                const string path,
                const string name,
                const mat &pos,
                const mat &vel,
                const vec &type)
{
    //Create an object for writing. The destructor will close the file.
    lammpswriter writer(nParticleProperties, frameNumber, nParticles, path, name);

    for (uint i = 0; i < nParticles; ++i)
    {
        writer << type[i%2]
               << pos(i, 0)
               << pos(i, 1)
               << pos(i, 2)
               << vel(i, 0)
               << vel(i, 1)
               << vel(i, 2);
    }

    writer.finalize();

}
