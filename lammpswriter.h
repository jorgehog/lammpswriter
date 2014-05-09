#pragma once

#include <fstream>
#include <sstream>
#include <exception>

#ifdef USE_MPI
#include <mpi.h>
#endif

using std::ofstream;
using std::string;
using std::cout;
using std::endl;

class lammpswriter
{

public:

    lammpswriter()
    {

    }

    lammpswriter(const uint frameNumber, const uint nParticles) :
        m_nParticles(nParticles),
        m_frameNumber(frameNumber)
    {
        initializeFile();

        dumpHeader();

    }

    ~lammpswriter()
    {
        finalize();
    }

    void initializeNewFile(const uint frameNumber, const uint nParticles)
    {
        m_frameNumber = frameNumber;

        m_nParticles = nParticles;

#ifndef NDEBUG
        _checkIfFileOpen();
#endif

        initializeFile();

        dumpHeader();
    }

    void finalize()
    {

#ifndef NDEBUG
        _checkParticlePropertySize();
#endif
        _safeGuardCounter = 0;
        m_file.close();

    }

    static void setSystemSize(const double systemSizeX,
                              const double systemSizeY,
                              const double systemSizeZ,
                              const double systemSizeX_start = 0,
                              const double systemSizeY_start = 0,
                              const double systemSizeZ_start = 0)
    {
        m_systemSizeX = systemSizeX;
        m_systemSizeY = systemSizeY;
        m_systemSizeZ = systemSizeZ;

        m_systemSizeX_start = systemSizeX_start;
        m_systemSizeY_start = systemSizeY_start;
        m_systemSizeZ_start = systemSizeZ_start;
    }


    static void setShear(const double xShear = 0,
                         const double yShear = 0,
                         const double zShear = 0)
    {
        m_xShear = xShear;
        m_yShear = yShear;
        m_zShear = zShear;
    }


    static void setNumberOfParticleProperties(const uint nParticleProperties)
    {
        m_nParticleProperties = nParticleProperties;
    }


    static void setFileNamePrefix(const string prefix)
    {
        m_prefix = prefix;
    }

    static void setFilePath(const string path)
    {
        m_path = path + "/";
    }

    static void setMPIRank(const int rank, const int masterRank)
    {
        m_MPIRank = rank;

        m_isMPIMaster = (rank == masterRank);
    }


    template<typename T, typename ...Args>
    void write(const T &val, const Args&... args)
    {
        write<T>(val);
        write(args...);
    }

    template<typename T>
    void write(const T &val)
    {
#ifndef NDEBUG
        _checkIfFileClosed();
        _safeGuardCounter++;
#endif
        m_file.write(reinterpret_cast<const char*>(&val), sizeof(val));
    }

    template<typename T>
    lammpswriter &operator << (const T &val)
    {
        write(static_cast<double>(val));
        return *this;
    }


private:

    static uint _safeGuardCounter;

    static double m_systemSizeX_start;
    static double m_systemSizeY_start;
    static double m_systemSizeZ_start;

    static double m_systemSizeX;
    static double m_systemSizeY;
    static double m_systemSizeZ;

    static double m_xShear;
    static double m_yShear;
    static double m_zShear;

    static uint m_nParticleProperties;

    static string m_prefix;
    static string m_path;

    static int m_MPIRank;

    static bool m_isMPIMaster;

    ofstream m_file;

    uint m_nParticles;

    uint m_frameNumber;


    void initializeFile()
    {
        std::stringstream s;

        s << m_path << m_prefix << m_frameNumber;

#ifdef USE_MPI
        s << "_" << m_MPIRank;
#endif
        s  << ".lmp";

        m_file.open(s.str().c_str());

    }

    void dumpHeader()
    {

        if (!m_isMPIMaster)
        {
            return;
        }

#ifndef NDEBUG
        _checkSystemSize();
#endif

        const uint chunkLength = m_nParticles*m_nParticleProperties;
        const uint nChunks = 1;

        write(m_frameNumber,
              m_nParticles,
              m_systemSizeX_start,
              m_systemSizeX,
              m_systemSizeY_start,
              m_systemSizeY,
              m_systemSizeZ_start,
              m_systemSizeZ,
              m_xShear,
              m_yShear,
              m_zShear,
              m_nParticleProperties,
              nChunks,
              chunkLength);

        _safeGuardCounter = 0;
    }


    void _checkSystemSize()
    {
        if (m_systemSizeX_start >= m_systemSizeX ||
                m_systemSizeY_start >= m_systemSizeY ||
                m_systemSizeZ_start >= m_systemSizeZ)
        {
            throw std::runtime_error("Inconsistent system sizes.");
        }
    }

    void _checkIfFileOpen()
    {
        if (m_file.is_open())
        {
            throw std::runtime_error("lammps file is already open. (Forgot to call finalize()?)");
        }
    }

    void _checkIfFileClosed()
    {
        if (!m_file.is_open())
        {
            throw std::runtime_error("lammps file is not open. (Forgot to call initializeNewFile()?)");
        }
    }

    void _checkParticlePropertySize()
    {
        if (_safeGuardCounter == 0)
        {
            return;
        }

        const uint nParticlePropertiesSaved = _safeGuardCounter/m_nParticles;
        const double _check = _safeGuardCounter/(double)m_nParticles;

        if (nParticlePropertiesSaved != _check)
        {
            throw std::runtime_error("Uneven number of particle properties saved.");
        }

        if (nParticlePropertiesSaved != m_nParticleProperties)
        {
            throw std::runtime_error("Saved number of particle properties does not match the specified number");
        }
    }

};

uint lammpswriter::_safeGuardCounter = 0;

double lammpswriter::m_systemSizeX_start;
double lammpswriter::m_systemSizeY_start;
double lammpswriter::m_systemSizeZ_start;

double lammpswriter::m_systemSizeX;
double lammpswriter::m_systemSizeY;
double lammpswriter::m_systemSizeZ;

double lammpswriter::m_xShear = 0;
double lammpswriter::m_yShear = 0;
double lammpswriter::m_zShear = 0;

uint lammpswriter::m_nParticleProperties;

string lammpswriter::m_prefix = "lammpsfile";
string lammpswriter::m_path = "";

int lammpswriter::m_MPIRank;
bool lammpswriter::m_isMPIMaster = true;
