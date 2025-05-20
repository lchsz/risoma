#include "utils.h"

using namespace std;

string trim(const string &str)
{
  size_t start = str.find_first_not_of(WHITESPACE);
  string lstr = (start == string::npos) ? "" : str.substr(start);
  size_t end = lstr.find_last_not_of(WHITESPACE);
  return (end == string::npos) ? "" : lstr.substr(0, end + 1);
}

int pack_ch(char ch)
{
  switch (ch)
  {
  case 'A':
    return BASE_A;
  case 'C':
    return BASE_C;
  case 'G':
    return BASE_G;
  case 'T':
    return BASE_T;
  }

  return 0;
}

vector<uint64_t> pack_read(const string &read_seq, int seed_len,
                           uint64_t mask)
{
  vector<uint64_t> out;

  int k_len = seed_len - 1;
  out.resize(read_seq.length() - k_len, 0);
  uint64_t pack = 0;
  for (int i = 0; i < read_seq.length(); i++)
  {
    uint64_t tmp = pack_ch(read_seq[i]);
    pack = (pack << 4) | tmp;
    if (i >= k_len)
    {
      pack &= mask;
      out[i - k_len] = pack;
    }
  }

  return out;
}

uint64_t pack_seed(const string &seed_seq)
{
  uint64_t pack = 0;

  for (int i = 0; i < seed_seq.length(); i++)
  {
    uint64_t tmp = pack_ch(seed_seq[i]);
    pack = (pack << 4) | tmp;
  }

  return pack;
}

int bit_count(uint64_t n)
{
  n = n - ((n >> 1) & 0x5555555555555555);
  n = (n & 0x3333333333333333) + ((n >> 2) & 0x3333333333333333);
  n = (n + (n >> 4)) & 0x0F0F0F0F0F0F0F0F;
  n = n + (n >> 8);
  n = n + (n >> 16);
  n = n + (n >> 32);
  return static_cast<int>(n & 0x7F);
}

// [[Rcpp::export]]
std::string generate_cigar(const std::string &seq1, const std::string &seq2)
{
  string cigar;
  int len = seq1.length();
  int count = 0;
  char last_op = '\0';

  for (int i = 0; i < len; i++)
  {
    char op;

    if (seq1[i] == '-' && seq2[i] == '-')
      continue;

    if (seq1[i] == '-' && seq2[i] != '-')
    {
      op = 'D';
    }
    else if (seq1[i] != '-' && seq2[i] == '-')
    {
      op = 'I';
    }
    else if (seq1[i] == seq2[i])
    {
      op = 'M';
    }
    else
    {
      op = 'X';
    }

    if (op == last_op)
    {
      ++count;
    }
    else
    {
      if (count > 0)
      {
        cigar += to_string(count) + last_op;
      }
      count = 1;
      last_op = op;
    }
  }

  if (count > 0)
  {
    cigar += to_string(count) + last_op;
  }

  return cigar;
}

int calc_snp_5p(const string &s1, const string &s2)
{
  int diff = 0;
  int i = s1.length() - 1;
  int j = s2.length() - 1;
  while (i >= 0 && j >= 0)
  {
    if (s1[i] != s2[j])
      ++diff;

    i--;
    j--;
  }

  return diff;
}

int calc_snp_3p(const string &s1, const string &s2)
{
  int diff = 0;
  int m = s1.length();
  int n = s2.length();
  int i = 0;
  int j = 0;
  while (i < m && j < n)
  {
    if (s1[i] != s2[j])
      ++diff;

    i++;
    j++;
  }

  return diff;
}

int64_t calc_mask(int length)
{
  int64_t mask = 0xFFFFFFFFFFFFFFFFL;

  if (length < 16)
  {
    mask <<= (16 - length) * 4;
  }

  return mask;
}

bool contains_n(const string &seq)
{
  return seq.find('N') != string::npos;
}

unordered_map<string, double>
calc_tpm(const unordered_map<string, int> &read_counts)
{
  unordered_map<string, double> rpk;
  double total_rrp = 0.0;

  for (const auto &entry : read_counts)
  {
    const string &gene = entry.first;
    int count = entry.second;
    int len = gene.length();
    double rpk_value = count / (len / 1000.0);
    rpk[gene] = rpk_value;
    total_rrp += rpk_value;
  }

  unordered_map<string, double> tpm;
  double scaling_factor = total_rrp / 1000000.0;

  for (const auto &entry : rpk)
  {
    const string &gene = entry.first;
    tpm[gene] = entry.second / scaling_factor;
  }

  return tpm;
}
