#include "utils.h"

string trim(const string &str)
{
  size_t start = str.find_first_not_of(WHITESPACE);
  string lstr = (start == string::npos) ? "" : str.substr(start);
  size_t end = lstr.find_last_not_of(WHITESPACE);
  return (end == string::npos) ? "" : lstr.substr(0, end + 1);
}