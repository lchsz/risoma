#include <Rcpp.h>
#include <stdio.h>
#include <zlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include "utils.h"

using namespace Rcpp;
using namespace std;

unordered_map<string, double> filter_reads(unordered_map<string, double>& read2tpm) {
  unordered_map<string, double> filtered_read2tpm;

  unordered_map<string, double>::const_iterator it;
  for (it = read2tpm.begin(); it != read2tpm.end(); it++)
  {
    const string& read_seq = it->first;
    double tpm = it->second;
    int read_len = read_seq.length();
    if (read_len >= 18 && read_len <= 26 && !contains_n(read_seq))
    {
      filtered_read2tpm[read_seq] = tpm;
    }
  }

  return filtered_read2tpm;
}


//' @title Load reads in a FASTQ File
//' @description This function reads a FASTQ file and marks duplicate reads
//'   based on their sequences.
//' @param fq_file Path to the FASTQ file.
//' @return A DataFrame containing two columns:
//'   - `read_seq`: The sequence of the read.
//'   - `tpm`: TPM of the read in the file.
//' @examples
//' \dontrun{
//'   result <- load_reads("example.fq")
//' }
//' @export
// [[Rcpp::export]]
DataFrame load_reads(std::string fq_file)
{
  ifstream fin;
  fin.open(fq_file.c_str());
  if (!fin.is_open())
  {
    stop("Could not open %s\n", fq_file.c_str());
  }

  unordered_map<string, int> read2num;

  string line;
  while (getline(fin, line))
  {
    if (line.substr(0, 1) == "@")
    {
      getline(fin, line);
      string read = trim(line);
      unordered_map<string, int>::iterator it = read2num.find(read);
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

  unordered_map<string, double> read2tpm = calc_tpm(read2num);
  unordered_map<string, double> filtered_read2tpm = filter_reads(read2tpm);

  vector<string> read_seqs;
  vector<double> tpms;

  unordered_map<string, double>::const_iterator it;
  for (it = filtered_read2tpm.begin(); it != filtered_read2tpm.end(); it++)
  {
      read_seqs.push_back(it->first);
      tpms.push_back(it->second);
  }

  DataFrame result = DataFrame::create(_["read_seq"] = read_seqs,
                                       _["tpm"] = tpms);

  return result;
}


//' @title Load reads in a GZIP-compressed FASTQ File
//' @description This function reads a GZIP-compressed FASTQ file and marks
//'   duplicate reads based on their sequences.
//' @param fq_file Path to the GZIP-compressed FASTQ file.
//' @return A DataFrame containing two columns:
//'   - `read_seq`: The sequence of the read.
//'   - `tpm`: TPM of the read in the file.
//' @examples
//' \dontrun{
//'   result <- load_read_gz("example.fq.gz")
//' }
//' @export
// [[Rcpp::export]]
DataFrame load_reads_gz(std::string fq_file)
{
  unordered_map<string, int> read2num;

  gzFile gz_file = gzopen(fq_file.c_str(), "r");
  if (NULL == gz_file)
  {
    printf("Could not open %s\n", fq_file.c_str());
    exit(EXIT_FAILURE);
  }

  char line[256];
  memset(line, 0, sizeof(line));
  while (gzgets(gz_file, line, sizeof(line)))
  {
    if (line[0] == '@')
    {
      gzgets(gz_file, line, sizeof(line));
      string read = trim(string(line));
      unordered_map<string, int>::iterator it = read2num.find(read);
      if (it != read2num.end())
      {
        it->second = it->second + 1;
      }
      else
      {
        read2num[read] = 1;
      }

      gzgets(gz_file, line, sizeof(line));
      gzgets(gz_file, line, sizeof(line));
    }
  }

  gzclose(gz_file);

  unordered_map<string, double> read2tpm = calc_tpm(read2num);
  unordered_map<string, double> filtered_read2tpm = filter_reads(read2tpm);

  vector<string> read_seqs;
  vector<double> tpms;

  unordered_map<string, double>::const_iterator it;
  for (it = filtered_read2tpm.begin(); it != filtered_read2tpm.end(); it++)
  {
      read_seqs.push_back(it->first);
      tpms.push_back(it->second);
  }

  DataFrame result = DataFrame::create(_["read_seq"] = read_seqs,
                                       _["tpm"] = tpms);

  return result;
}


