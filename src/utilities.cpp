#include "utilities.hpp"

#include <cstring>
#include <fstream>


int argcheck(int argc, char **argv, const char *arg)
{
  for(int i = 1; i < argc; ++i)
  {
    // strcmp returns 0 if the strings are equal
    if(strcmp(argv[i], arg) == 0)
      return(i);
  }
  return 0;
}



void read_binary(const std::string &filename, int n_values, double *values)
{
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in)
    throw std::runtime_error("File '" + filename + "' can't be opened");

  in.seekg(0, in.end); // jump to the end of the file
  int length = in.tellg(); // total length of the file in bytes
  int size_value = length / n_values; // size (in bytes) of one value

  if (length % n_values != 0)
    throw std::runtime_error("The number of bytes in the file '" + filename +
                             "' is not divisible by the number of elements " +
                             d2s(n_values));

  in.seekg(0, in.beg); // jump to the beginning of the file

  if (size_value == sizeof(double))
  {
    // read all at once
    in.read(reinterpret_cast<char*>(values), n_values*size_value);

    if(n_values != static_cast<int>(in.gcount()))
      throw std::runtime_error("The number of successfully read elements is "
                               "different from the expected one");
  }
  else if (size_value == sizeof(float))
  {
    float val;
    for (int i = 0; i < n_values; ++i)  // read element-by-element
    {
      // read a 'float' value and convert it to a 'double' value
      in.read(reinterpret_cast<char*>(&val), size_value);
      values[i] = val;
    }
  }
  else
    throw std::runtime_error("Unknown size of an element (" + d2s(size_value) +
                             ") in bytes. Expected one is either sizeof(float) "
                             "= " + d2s(sizeof(float)) + ", or sizeof(double) "
                             "= " + d2s(sizeof(double)));
}



void write_binary(const std::string &filename, int n_values,
                  const double *values)
{
  std::ofstream out(filename.c_str(), std::ios::binary);
  if (!out)
    throw std::runtime_error("File '" + filename + "' can't be opened");

  for (int i = 0; i < n_values; ++i)
  {
    float val = values[i];
    out.write(reinterpret_cast<char*>(&val), sizeof(float));
  }
}



std::string file_name(const std::string &path)
{
  if (path == "") return path;
  // extract a filename
  const std::string fname = path.substr(path.find_last_of('/') + 1);
  return fname;
}



std::string file_stem(const std::string &path)
{
  if (path == "") return path;
  // get a file name from the path
  const std::string fname = file_name(path);
  // extract a stem and return it
  return fname.substr(0, fname.find_last_of('.'));
}

