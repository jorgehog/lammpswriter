#pragma once

#include <sys/types.h>
#include <vector>

using std::vector;

namespace _lammpswriterDataHandler
{

class Base
{
public:

    Base() {}

    virtual ~Base()
    {
        m_myValues.clear();
    }

    virtual void push(const double value, const uint index) = 0;
    virtual void initialize() = 0;

    uint size() const
    {
        return m_myValues.size();
    }

    const vector<double> &myValues() const
    {
        return m_myValues;
    }

protected:

    vector<double> m_myValues;

};


class UndeterminedSize : public Base
{
public:

    UndeterminedSize() : Base() {}


    void push(const double value, const uint index)
    {
        (void) index;
        m_myValues.push_back(value);
    }

    void initialize()
    {

    }

};


class DeterminedSize : public Base
{
public:

    DeterminedSize(const uint N) : Base(), m_N(N) {}

    void push(const double value, const uint index)
    {
        m_myValues[index] = value;
    }

    void initialize()
    {
        m_myValues.resize(m_N);
    }
    void finalize()
    {
    }

private:

    const uint m_N;

    void _checkBounds(const uint &index)
    {
#ifndef NDEBUG
        if (index >= m_N)
        {
            throw std::runtime_error("Index out of bounds: More properties stored than specified.");
        }
#else
        (void) index;
#endif
    }

};

}
