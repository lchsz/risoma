#include <Rcpp.h>
#include <iostream>
#include <fstream>

using namespace Rcpp;
using namespace std;

const string WHITESPACE = " \n\r\t\f\v";

string trim(const string &str)
{
  size_t start = str.find_first_not_of(WHITESPACE);
  string lstr = (start == string::npos) ? "" : str.substr(start);
  size_t end = lstr.find_last_not_of(WHITESPACE);
  return (end == string::npos) ? "" : lstr.substr(0, end + 1);
}

//' Mark duplicates in the FASTQ file
//'
//' To speed up the computation, we first extract sequences of each read and
//' then mark and count duplicates.
//'
//' @param fqFile FASTQ file of a sample
//' @param minReadNum the minimum read count for an isoform
//' @return a data.frame with column 'read_seq' and 'read_num'
// [[Rcpp::export]]
DataFrame markDuplicates(std::string fqFile, int minReadNum)
{
  map<string, int> read2num;

  ifstream fin;
  fin.open(fqFile.c_str());
  if (!fin.is_open())
  {
    stop("Could not open %s\n", fqFile.c_str());
  }

  string line;
  while (getline(fin, line))
  {
    if (line.substr(0, 1) == "@")
    {
      getline(fin, line);
      string read = trim(line);
      map<string, int>::iterator it = read2num.find(read);
      if (it != read2num.end())
      {
        it->second = it->second + 1;
      }
      else
      {
        read2num[read] = 1;
      }
      getline(fin, line);
      getline(fin, line);
    }
  }
  fin.close();

  vector<string> readSeqs;
  vector<int> readNums;

  map<string, int>::const_iterator it;
  for (it = read2num.begin(); it != read2num.end(); it++)
  {
    if (it->second > minReadNum)
    {
      readSeqs.push_back(it->first);
      readNums.push_back(it->second);
    }
  }

  DataFrame result = DataFrame::create(_["read_seq"] = readSeqs,
                                       _["read_num"] = readNums);

  return result;
}
