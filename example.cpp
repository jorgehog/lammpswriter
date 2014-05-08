#include <iostream>
#include <sys/types.h>
#include <armadillo>

#include "lammpswriter/lammpswriter.h"

using namespace arma;
using namespace std;

void dumpLammps(const uint frameNumber, const uint nParticles, const mat &pos, const mat &vel, const vec &type);

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

    //Set the system size once and for all.
    lammpswriter::setSystemSize(1, 1, 1);

    //Set the number of properties to store for each particle.
    lammpswriter::setNumberOfParticleProperties(7);

    //Set a path to save the files. Defaults to current working directory.
    lammpswriter::setFilePath("/tmp");

    //Set the name of the files (particlePositions0.lmp ...).
    //Defaults to 'lammpsfile'
    lammpswriter::setFileNamePrefix("particlePositions");


    uint nParticles = 10;

    //We store positions, velocities, and a type
    mat pos(nParticles, 3);
    mat vel(nParticles, 3);
    vec type = {1.0, 2.0};


    //Create an object for writing (used only for alternative one)
    lammpswriter writer;


    uint nCycles = 10;

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
            writer.write(type[i%2],
                         pos(i, 0),
                         pos(i, 1),
                         pos(i, 2),
                         vel(i, 0),
                         vel(i, 1),
                         vel(i, 2));
        }

        //closes the file.
        writer.finalize();
    }

    //Second alternative usage
    for (uint c = nCycles; c < 2*nCycles; ++c)
    {
        pos.randu();
        vel.randn();

        dumpLammps(c, nParticles, pos, vel, type);
    }

    return 0;
}

void dumpLammps(const uint frameNumber, const uint nParticles, const mat &pos, const mat &vel, const vec &type)
{
    //Create an object for writing. The destructor will close the file.
    lammpswriter writer(frameNumber, nParticles);

    for (uint i = 0; i < nParticles; ++i)
    {
        writer.write(type[i%2],
                     pos(i, 0),
                     pos(i, 1),
                     pos(i, 2),
                     vel(i, 0),
                     vel(i, 1),
                     vel(i, 2));
    }
}
