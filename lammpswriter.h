#pragma once

#include <fstream>
#include <sstream>
#include <exception>
#include <vector>
#include <numeric>

#ifdef LAMMPSWRITER_USE_MPI
#include <mpi.h>
#endif

using std::ifstream;
using std::ofstream;
using std::string;
using std::cout;
using std::endl;
using std::vector;

class lammpswriter
{

public:

    lammpswriter(const uint nParticleProperties,
                 const string prefix = "lammpsfile",
                 const string path = ""):
        m_systemSizeX_start(0),
        m_systemSizeY_start(0),
        m_systemSizeZ_start(0),

        m_systemSizeX(0),
        m_systemSizeY(0),
        m_systemSizeZ(0),

        m_xShear(0),
        m_yShear(0),
        m_zShear(0),

        m_valueCounter(0),
        m_nParticleProperties(nParticleProperties),
        m_path(path.empty() ? path : path + "/"),
        m_prefix(prefix)
    {

    }

    ~lammpswriter()
    {
        _checkIfFileOpen();

        m_myValues.clear();
        m_allValues.clear();

    }

    void initializeNewFile(const uint frameNumber, const uint nParticles)
    {
        m_fileState = OUT;

        m_frameNumber = frameNumber;

        m_nParticles = nParticles;

        _checkIfFileOpen();

        _initialize();

    }

    void finalize()
    {

        _checkParticlePropertySize();

        if (m_fileState == OUT)
        {
            _dumpFile();
        }
        else
        {
            m_inFile.close();
        }

        m_valueCounter = 0;

    }

    void loadFile(const uint framenNumber)
    {
        m_fileState = IN;

        m_frameNumber = framenNumber;

        m_inFile.open(_fileName(), std::ios::binary);

        _checkIfFileExists();

        uint frameNumber, nChunks, chunkLength, nParticleProperties;

        read(frameNumber,
             m_totalParticles,
             m_systemSizeX_start,
             m_systemSizeX,
             m_systemSizeY_start,
             m_systemSizeY,
             m_systemSizeZ_start,
             m_systemSizeZ,
             m_xShear,
             m_yShear,
             m_zShear,
             nParticleProperties,
             nChunks,
             chunkLength);

        m_nParticles = m_totalParticles;

        if (frameNumber != m_frameNumber)
        {
            throw std::logic_error("Mismatch in specified and loaded frame.");
        }

        if (nParticleProperties != m_nParticleProperties)
        {
            throw std::logic_error("Mismatch in specified and loaded number of columns.");
        }

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
        _checkIfFileClosed();

        m_file.write(reinterpret_cast<const char*>(&val), sizeof(val));
    }

    template<typename T, typename ...Args>
    void read(T &val, Args&... args)
    {
        read<T>(val);
        read(args...);
    }

    template<typename T>
    void read(T &val)
    {
        _checkIfFileClosed();

        m_inFile.read(reinterpret_cast<char*>(&val), sizeof(T));
    }

    template<typename T>
    lammpswriter &operator << (const T &val)
    {
        m_myValues[m_valueCounter++] = static_cast<double>(val);
        return *this;
    }

    template<typename T>
    lammpswriter &operator >> (T &val)
    {
        double valueAsDouble = static_cast<double>(val);

        m_inFile.read(reinterpret_cast<char*>(&valueAsDouble), sizeof(double));

        m_valueCounter++;

        val = static_cast<T>(valueAsDouble);

        return *this;
    }

    lammpswriter &operator >> (double &val)
    {
        m_inFile.read(reinterpret_cast<char*>(&val), sizeof(double));

        m_valueCounter++;

        return *this;
    }

    void setSystemSize(const double systemSizeX,
                       const double systemSizeY,
                       const double systemSizeZ,
                       const double systemSizeX_start = 0,
                       const double systemSizeY_start = 0,
                       const double systemSizeZ_start = 0)
    {
        _checkIfFileOpen();

        m_systemSizeX = systemSizeX;
        m_systemSizeY = systemSizeY;
        m_systemSizeZ = systemSizeZ;

        m_systemSizeX_start = systemSizeX_start;
        m_systemSizeY_start = systemSizeY_start;
        m_systemSizeZ_start = systemSizeZ_start;
    }


    void setShear(const double xShear = 0,
                  const double yShear = 0,
                  const double zShear = 0)
    {
        _checkIfFileOpen();

        m_xShear = xShear;
        m_yShear = yShear;
        m_zShear = zShear;
    }

    static void setMPIRank(const int rank, const int nNodes, const int masterRank = 0)
    {
        m_MPI_master = masterRank;

        m_MPI_nNodes = nNodes;

        _checkMPI();

        m_isMPIMaster = (rank == masterRank);

        if (m_isMPIMaster)
        {
            m_nParticlesList.resize(nNodes);
        }

    }

    void setPath(const string path)
    {
        _checkIfFileOpen();

        m_path = path.empty() ? path : path + "/";
    }

    void setPrefix(const string prefix)
    {
        _checkIfFileOpen();

        m_prefix = prefix;
    }

    void setNParticleProperties(const uint nParticleProperties)
    {
        _checkIfFileOpen();

        m_nParticleProperties = nParticleProperties;
    }


    const double &systemSizeX_start() const
    {
        return m_systemSizeX_start;
    }

    const double &systemSizeY_start() const
    {
        return m_systemSizeY_start;
    }

    const double &systemSizeZ_start() const
    {
        return m_systemSizeZ_start;
    }


    const double &systemSizeX() const
    {
        return m_systemSizeX;
    }

    const double &systemSizeY() const
    {
        return m_systemSizeX;
    }

    const double &systemSizeZ() const
    {
        return m_systemSizeX;
    }


    const double &xShear() const
    {
        return m_xShear;
    }

    const double &yShear() const
    {
        return m_yShear;
    }

    const double &zShear() const
    {
        return m_zShear;
    }

    const uint &nParticleProperties() const
    {
        return m_nParticleProperties;
    }

    const string &path() const
    {
        return m_path;
    }

    const string &prefix() const
    {
        return m_prefix;
    }

    const uint &nParticles() const
    {
        return m_nParticles;
    }
    const uint &totalParticles() const
    {
        return m_totalParticles;
    }


private:

    double m_systemSizeX_start;
    double m_systemSizeY_start;
    double m_systemSizeZ_start;

    double m_systemSizeX;
    double m_systemSizeY;
    double m_systemSizeZ;

    double m_xShear;
    double m_yShear;
    double m_zShear;


    static int m_MPI_master;
    static bool m_isMPIMaster;
    static int m_MPI_nNodes;
    static vector<int> m_nParticlesList;

    enum fileState
    {
        IN,
        OUT
    } m_fileState;


    uint m_valueCounter;

    uint m_nParticleProperties;

    string m_path;
    string m_prefix;

    vector<double> m_myValues;
    vector<double> m_allValues;

    ofstream m_file;
    ifstream m_inFile;

    uint m_nParticles;
    uint m_totalParticles;

    uint m_frameNumber;

    const char *_fileName() const
    {
        std::stringstream s;

        s << m_path << m_prefix << m_frameNumber << ".lmp";

        return s.str().c_str();
    }

    void _initialize()
    {

        m_myValues.resize(m_nParticles*m_nParticleProperties);

        if (!m_isMPIMaster)
        {
            return;
        }

        m_file.open(_fileName(), std::ios::binary);

        _checkIfFileExists();

        _dumpHeader();

    }

    void _dumpHeader()
    {

        _getTotalNumberOfParticles();

        if (!m_isMPIMaster)
        {
            return;
        }

        _checkSystemSize();

        const uint chunkLength = m_totalParticles*m_nParticleProperties;
        const uint nChunks = 1;

        write(m_frameNumber,
              m_totalParticles,
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
#ifndef NDEBUG
        if (m_file.is_open())
        {
            throw std::runtime_error("lammps file is already open. (Forgot to call finalize()?)");
        }
#endif
    }

    void _checkIfFileClosed()
    {
#ifndef NDEBUG
        if (m_fileState == OUT)
        {
            if (!m_file.is_open())
            {
                throw std::runtime_error("lammps file is not open. (Forgot to call initializeNewFile()?)");
            }
        }
        else
        {
            if (!m_inFile.is_open())
            {
                throw std::runtime_error("lammps file is not open. (Forgot to call initializeNewFile()?)");
            }
        }
#endif
    }


    void _checkIfFileExists()
    {
#ifndef NDEBUG
        string error = "lammps file could not be opened. Bad file? " + string(_fileName());

        if (m_fileState == OUT)
        {
            if (!m_file.good())
            {
                throw std::runtime_error(error);
            }
        }
        else
        {
            if (!m_inFile.good())
            {
                throw std::runtime_error(error);
            }
        }
#endif
    }

    void _checkParticlePropertySize()
    {
#ifndef NDEBUG
        if (m_valueCounter == 0)
        {
            return;
        }

        const uint nParticlePropertiesSaved = m_valueCounter/m_nParticles;
        const double _check = m_valueCounter/(double)m_nParticles;

        if (nParticlePropertiesSaved != _check)
        {
            throw std::runtime_error("Uneven number of particle properties saved.");
        }

        if (nParticlePropertiesSaved != m_nParticleProperties)
        {
            throw std::runtime_error("Saved number of particle properties does not match the specified number");
        }
#endif
    }

    static void _checkMPI()
    {
#ifndef NDEBUG
#ifdef LAMMPSWRITER_LOW_MEMORY
        if (m_MPI_master != 0)
        {
            throw std::runtime_error("For low mem mpi to work, master rank must be rank zero.");
        }
#endif
#endif
    }

    void _getTotalNumberOfParticles()
    {
#ifdef LAMMPSWRITER_USE_MPI

        MPI_Gather(&m_nParticles, 1, MPI_INT, &m_nParticlesList.front(), 1, MPI_INT, m_MPI_master, MPI_COMM_WORLD);

        if (m_isMPIMaster)
        {
            m_totalParticles = std::accumulate(m_nParticlesList.begin(), m_nParticlesList.end(), 0);
        }

#else
        m_totalParticles = m_nParticles;
#endif
    }

    void _dumpFile()
    {
#ifdef LAMMPSWRITER_USE_MPI
        if (m_MPI_nNodes != 1)
        {

            int *displacements = NULL;
            int *recvCounts = NULL;

            if (m_isMPIMaster)
            {

                displacements = new int[m_MPI_nNodes];
                recvCounts = new int[m_MPI_nNodes];

                int currentDisplacement = 0;
                for (int i = 0; i < m_MPI_nNodes; ++i)
                {
                    displacements[i] = currentDisplacement;

                    recvCounts[i] = m_nParticlesList[i]*m_nParticleProperties;

                    currentDisplacement += recvCounts[i];
                }
#ifndef LAMMPSWRITER_LOW_MEMORY
                m_allValues.resize(m_totalParticles*m_nParticleProperties);
#endif
            }

#ifndef LAMMPSWRITER_LOW_MEMORY
            MPI_Gatherv(&m_myValues.front(),
                        m_myValues.size(),
                        MPI_DOUBLE,
                        &m_allValues.front(),
                        recvCounts,
                        displacements,
                        MPI_DOUBLE,
                        m_MPI_master,
                        MPI_COMM_WORLD);
#else
            if (m_isMPIMaster)
            {
                _simpleDump();

                for (uint node = 1; node < m_MPI_nNodes; node++)
                {
                    m_myValues.resize(recvCounts[node]);
                    MPI_Recv(&m_myValues.front(), m_myValues.size(), MPI_DOUBLE, node, 0, MPI_COMM_WORLD, NULL);

                    m_nParticles = m_nParticlesList[node];
                    _simpleDump();
                }
            }

            else
            {
                MPI_Send(&m_myValues.front(), m_myValues.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
#endif


            if (m_isMPIMaster)
            {
#ifndef LAMMPSWRITER_LOW_MEMORY
                m_file.write(reinterpret_cast<const char*>(&m_allValues.front()), m_totalParticles*m_nParticleProperties*sizeof(double));
#endif

                delete [] recvCounts;
                delete [] displacements;
            }
        }

        else
        {
            _simpleDump();
        }

#else
        _simpleDump();
#endif
        m_file.close();
    }

    void _simpleDump()
    {
        m_file.write(reinterpret_cast<const char*>(&m_myValues.front()), m_nParticles*m_nParticleProperties*sizeof(double));
    }

};

int  lammpswriter::m_MPI_master = 0;
bool lammpswriter::m_isMPIMaster = true;
int  lammpswriter::m_MPI_nNodes;
vector<int>    lammpswriter::m_nParticlesList;
