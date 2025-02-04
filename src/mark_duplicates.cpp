#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "utils.h"

using namespace Rcpp;
using namespace std;

//' Mark duplicates in the FASTQ file
//'
//' To speed up the computation, we first extract sequences of each read and
//' then mark and count duplicates.
//'
//' @param fqFile FASTQ file of a sample
//' @param minReadNum the minimum read count for an isoform
//' @return a data.frame with column 'read_seq' and 'read_num'
// [[Rcpp::export]]
DataFrame mark_duplicates(std::string fq_file, int min_read_num)
{
  map<string, int> read2num;

  ifstream fin;
  fin.open(fq_file.c_str());
  if (!fin.is_open())
  {
    stop("Could not open %s\n", fq_file.c_str());
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

  vector<string> read_seqs;
  vector<int> read_nums;

  map<string, int>::const_iterator it;
  for (it = read2num.begin(); it != read2num.end(); it++)
  {
    if (it->second > min_read_num)
    {
      read_seqs.push_back(it->first);
      read_nums.push_back(it->second);
    }
  }

  DataFrame result = DataFrame::create(_["read_seq"] = read_seqs,
                                       _["read_num"] = read_nums);

  return result;
}
