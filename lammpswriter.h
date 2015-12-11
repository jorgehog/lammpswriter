#pragma once

#include "datahandler.h"

#include <sys/types.h>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <numeric>
#include <cmath>

#ifdef LAMMPSWRITER_USE_MPI
#include <mpi.h>
#endif

using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;

class lammpswriter
{

public:

    lammpswriter(const uint nParticleProperties,
                 const string prefix = "lammpsfile",
                 const string path = ""):
        m_systemSizeSet(false),

        m_systemSizeX_start(0),
        m_systemSizeY_start(0),
        m_systemSizeZ_start(0),

        m_systemSizeX(0),
        m_systemSizeY(0),
        m_systemSizeZ(0),

        m_xShear(0),
        m_yShear(0),
        m_zShear(0),

        m_MPI_master(0),
        m_isMPIMaster(true),

        m_valueCounter(0),
        m_nParticleProperties(nParticleProperties),
        m_path(path.empty() ? path : path + "/"),
        m_prefix(prefix),

        m_dataHandler(nullptr)
    {

    }

    ~lammpswriter()
    {
        _checkIfFileOpen();

        m_allValues.clear();

    }

    void initializeNewFile(const uint frameNumber, const uint nParticles)
    {
       _initializeFrame(frameNumber, new _lammpswriterDataHandler::DeterminedSize(nParticles*m_nParticleProperties));
    }

    void initializeNewFile(const uint frameNumber)
    {
       _initializeFrame(frameNumber, new _lammpswriterDataHandler::UndeterminedSize());
    }

    void initializeNewFile()
    {
        initializeNewFile(0);
    }

    void finalize()
    {

        _checkParticlePropertySize();

        m_nParticles = m_valueCounter / m_nParticleProperties;

        if (m_fileState == OUT)
        {
            _dumpHeader();

            _dumpFile();

            delete m_dataHandler;
        }
        else
        {
            m_inFile.close();
        }

        m_valueCounter = 0;

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
        _checkIfFileClosed();

        m_dataHandler->push(static_cast<double>(val), m_valueCounter++);

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


    void loadFile(const uint framenNumber)
    {
        m_fileState = IN;

        m_frameNumber = framenNumber;

        m_inFile.open(_fileName().c_str(), std::ios::binary);

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

        if (frameNumber != m_frameNumber)
        {
            throw std::logic_error("Mismatch in specified and loaded frame.");
        }

        if (nParticleProperties != m_nParticleProperties)
        {
            throw std::logic_error("Mismatch in specified and loaded number of columns.");
        }

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

        m_systemSizeSet = true;
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

    void setMPIRank(const int rank, const int nNodes, const int masterRank = 0)
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
        return m_systemSizeY;
    }

    const double &systemSizeZ() const
    {
        return m_systemSizeZ;
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

    uint nParticles() const
    {
        return m_dataHandler->size()/m_nParticleProperties;
    }

    const uint &totalParticles() const
    {
        return m_totalParticles;
    }

    const uint &valueCounter() const
    {
        return m_valueCounter;
    }


private:

    bool m_systemSizeSet;

    double m_systemSizeX_start;
    double m_systemSizeY_start;
    double m_systemSizeZ_start;

    double m_systemSizeX;
    double m_systemSizeY;
    double m_systemSizeZ;

    double m_xShear;
    double m_yShear;
    double m_zShear;

    int m_MPI_master;
    bool m_isMPIMaster;
    int m_MPI_nNodes;
    vector<int> m_nParticlesList;

    enum fileState
    {
        IN,
        OUT
    } m_fileState;

    uint m_valueCounter;

    uint m_nParticleProperties;

    string m_path;
    string m_prefix;

    vector<double> m_allValues;

    _lammpswriterDataHandler::Base *m_dataHandler;

    ofstream m_file;
    ifstream m_inFile;

    uint m_nParticles;
    uint m_totalParticles;

    uint m_frameNumber;

    string _fileName() const
    {
        std::stringstream s;

        s << m_path << m_prefix << m_frameNumber << ".lmp";

        return s.str();
    }

    void _initializeFrame(const uint frameNumber, _lammpswriterDataHandler::Base *dataHandler)
    {
        m_fileState = OUT;

        m_frameNumber = frameNumber;

        m_dataHandler = dataHandler;

        _checkIfFileOpen();

        _initialize();
    }

    void _initialize()
    {

        m_dataHandler->initialize();

        if (!m_isMPIMaster)
        {
            return;
        }

        m_file.open(_fileName().c_str(), std::ios::binary);

        _checkIfFileExists();

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





    void _dumpFile()
    {
#ifdef LAMMPSWRITER_USE_MPI
        if (m_MPI_nNodes != 1)
        {

            int *displacements = nullptr;
            int *recvCounts = nullptr;

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
            MPI_Gatherv(&m_dataHandler->myValues().front(),
                        m_dataHandler->size(),
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

                for (int node = 1; node < m_MPI_nNodes; node++)
                {
                    m_dataHandler->resize(recvCounts[node]);
                    MPI_Recv(&m_dataHandler->myValues().front(), m_dataHandler->size(), MPI_DOUBLE, node, 0, MPI_COMM_WORLD, NULL);

                    m_nParticles = m_nParticlesList[node];
                    _simpleDump();
                }
            }

            else
            {
                MPI_Send(&m_dataHandler->myValues().front(), m_dataHandler->size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
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
        m_file.write(reinterpret_cast<const char*>(&m_dataHandler->myValues().front()), m_dataHandler->size()*sizeof(double));
    }

    void _checkSystemSize()
    {
#ifndef NDEBUG
        if (!m_systemSizeSet)
        {
            throw std::runtime_error("System size not set. Forgot to call setSystemSize() before initialize?");
        }

        else if (m_systemSizeX_start > m_systemSizeX ||
                 m_systemSizeY_start > m_systemSizeY ||
                 m_systemSizeZ_start > m_systemSizeZ)
        {
            throw std::runtime_error("Inconsistent system sizes.");
        }
#endif
    }

    void _checkIfFileOpen()
    {
#ifndef NDEBUG
        if (m_file.is_open())
        {
            throw std::runtime_error("lammps file is already open. Forgot to call finalize()?");
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
                throw std::runtime_error("lammps file is not open. Forgot to call initializeNewFile()?");
            }
        }
        else
        {
            if (!m_inFile.is_open())
            {
                throw std::runtime_error("lammps file is not open. Forgot to call initializeNewFile()?");
            }
        }
#endif
    }


    void _checkIfFileExists()
    {
#ifndef NDEBUG
        std::stringstream error;
        error << "lammps file could not be opened. Bad file? " << _fileName();

        if (m_fileState == OUT)
        {
            if (!m_file.good())
            {
                throw std::runtime_error(error.str());
            }
        }
        else
        {
            if (!m_inFile.good())
            {
                throw std::runtime_error(error.str());
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

        const double nParticlePropertiesSaved = m_valueCounter/(double)nParticles();
        const double _check = m_valueCounter/(double)nParticles();

        if (fabs(nParticlePropertiesSaved - (int)m_nParticleProperties) > 0.01)
        {
            std::stringstream s;
            s << "Saved number of particle properties does not match the specified number: " << nParticlePropertiesSaved << " is stored, and " << m_nParticleProperties << " is requested. ";

            throw std::runtime_error(s.str());
        }

        if (nParticlePropertiesSaved != _check)
        {
            std::stringstream s;
            s << "Uneven number of particle properties saved: " << nParticlePropertiesSaved << " is stored, and " << _check << " is requested. ";

            throw std::runtime_error(s.str());
        }

#endif
    }

    void _checkMPI()
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

        uint m_nParticles = nParticles();
        MPI_Gather(&m_nParticles, 1, MPI_INT, &m_nParticlesList.front(), 1, MPI_INT, m_MPI_master, MPI_COMM_WORLD);

        if (m_isMPIMaster)
        {
            m_totalParticles = std::accumulate(m_nParticlesList.begin(), m_nParticlesList.end(), 0);
        }

#else
        m_totalParticles = nParticles();
#endif
    }



};
