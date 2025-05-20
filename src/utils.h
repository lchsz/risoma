#include <string>
#include <cstdint>
#include <unordered_map>
#include <Rcpp.h>

const std::string WHITESPACE = " \n\r\t\f\v";

const int BASE_A = 0x1;
const int BASE_C = 0x4;
const int BASE_G = 0x8;
const int BASE_T = 0x2;

std::string trim(const std::string &str);

std::vector<uint64_t> pack_read(const std::string &read_seq, int seed_len,
                                uint64_t mask);

uint64_t pack_seed(const std::string &seq);

int bit_count(uint64_t n);

int calc_snp_5p(const std::string &s1, const std::string &s2);

int calc_snp_3p(const std::string &s1, const std::string &s2);

bool contains_n(const std::string &seq);

std::unordered_map<std::string, double>
calc_tpm(const std::unordered_map<std::string, int> &read_counts);
