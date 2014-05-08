#pragma once

#include <fstream>
#include <sstream>
#include <exception>

using std::ofstream;
using std::string;

class lammpswriter
{

public:

    lammpswriter(const uint frameNumber, const uint nParticles) :
        m_nParticles(nParticles),
        m_frameNumber(frameNumber)
    {
        m_chunkLength = m_nParticles*m_nParticleProperties;

        initializeFile();

        dumpHeader();
    }

    ~lammpswriter()
    {
        finalize();
    }


    static void setSystemSize(const uint systemSizeX,
                              const uint systemSizeY,
                              const uint systemSizeZ,
                              const uint systemSizeX_start = 0,
                              const uint systemSizeY_start = 0,
                              const uint systemSizeZ_start = 0)
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
        m_path = path;
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
        std::cout << "dumping " << val << " of size " << sizeof(val) << std::endl;
        m_file.write(reinterpret_cast<const char*>(&val), sizeof(val));
    }

private:

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

    static const uint m_nChunks;
    static  uint m_chunkLength;

    ofstream m_file;

    const uint &m_nParticles;

    const uint &m_frameNumber;


    void initializeFile()
    {
        std::stringstream s;
        s << m_path << "/" << m_prefix << m_frameNumber << ".lmp";

        m_file.open(s.str().c_str());
    }

    void dumpHeader()
    {
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
              m_nChunks,
              m_chunkLength);
    }

    void finalize()
    {
        m_file.close();
    }

};


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

const uint lammpswriter::m_nChunks = 1;
uint lammpswriter::m_chunkLength;
