#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <stdexcept>
#include <string>
#include <sstream>

/**
 * Convert the data of any type which has oveloaded operator '<<' to string
 * @return data in string format
 */
template <typename T>
inline std::string d2s(T data)
{
  std::ostringstream o;
  if (!(o << data))
    throw std::runtime_error("Bad conversion of data to string!");
  return o.str();
}

/**
 * @brief Test for presence of a command line argument
 * Compares a specified string to the array of standard Unix commnand line
 * arguments. If the string  matches the i'th argument, it returns i, otherwise
 * returns 0.
 * @author R. Gibson
 * @param argc - from command line, number of command line arguments
 * @param argv - pointer to array of command line arguments
 * @param arg - pointer to string to test for
 * @return i, if the string  matches the i'th argument, 0 otherwise.
 */
int argcheck(int argc, char **argv, const char *arg);

/**
 * Read an array of floating-point numbers in double precision from a binary
 * file.
 * @param filename - name of the binary file
 * @param n_values - number of values in the binary file
 * @param values - an array of values in double precision
 */
void read_binary(const std::string &filename, int n_values, double *values);

/**
 * Write a binary file with single precision floating point numbers
 * @param filename - name of the binary file
 * @param n_values - number of values to be written
 * @param values - values to be written
 */
void write_binary(const std::string &filename, int n_values,
                  const double *values);

/**
 * Get (extract) a file name from the given path.
 * @param path - a name of a file under interest including the path
 * @return a string representing a name of the file. For example:
 * @verbatim
   file_name("/home/user/file.dat") = "file.dat"
 * @endverbatim
 */
std::string file_name(const std::string &path);

/**
 * Extract a stem from a filename with a path.
 * @param path - a name of a file under interest including the path
 * @return a string which represents the name of the file without an extension -
 * only a stem of the file. For example:
 * @verbatim
   stem("/home/user/file.dat") = "file"
 * @endverbatim
 */
std::string file_stem(const std::string &path);

#endif // UTILITIES_HPP
