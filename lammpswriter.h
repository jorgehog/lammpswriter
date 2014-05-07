#pragma once

#include <fstream>

using std::ofstream;

namespace lammpswriter
{

template<typename T, typename ...Args>
static void write(ofstream &file, const T &val, const Args&... args)
{
    write<T>(file, val);
    write(file, args...);
}

template<typename T>
static void write(ofstream &file, const T &val)
{
    file.write(reinterpret_cast<const char*>(&val), sizeof(val));
}

}
