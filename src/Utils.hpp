#ifndef _UTILS_H
#define _UTILS_H

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>


using namespace std;

extern bool DEBUG;



/**
 * @brief BF Associates a Boltzmann factor to an energy contribution
 * @param dG Free-energy contribution
 * @return Boltzmann factor exp(-dG/RT) of contribution
 */
double BF(double dG);

/***************** Strings Utils *****************/
/**
 * @brief ltrim  Trims a string from the start
 * @param s String whose prefix must be trimmed
 * @return Prefix-trimmed string
 */
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

/**
 * @brief rtrim  Trims a string from the end
 * @param s String whose suffix must be trimmed
 * @return Suffix-trimmed string
 */
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

/**
 * @brief trim  Trims a string from both ends
 * @param s String to be trimmed
 * @return Trimmed string
 */
static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}

/**
 * @brief split Splits a string according to a given delimiter
 * @param s String to be split
 * @param delim Delimiter character
 * @return Vector of substrings split with respect to delimiter
 */
vector<string> split(const std::string &s, char delim);

/**
 * @brief pad Creates a string consisting of copies of a given character
 * @param c Character to copy
 * @param a Number of copies
 * @return String consisting of the concatenation of a copies of c
 */
string pad(const char * c, int a);

/**
 * @brief operator << Outputs to stream a given vector of integers
 * @param o Output stream
 * @param v Vector of integers
 * @return Output stream
 */
ostream& operator<<(ostream& o, const vector<int> & v);

/**************************************************/

/**
 * @brief setSubstract Basic set subtraction operation
 * @param indices First operand
 * @param parent Second operand
 * @return List of indices that are in indices but not in parent
 */
vector<int> setSubstract(vector<int> indices,vector<int> parent);

#endif 
