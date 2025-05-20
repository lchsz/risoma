#include <Rcpp.h>
#include "utils.h"
#include <bitset>

using namespace Rcpp;
using namespace std;

// @title Isoform Structure
// @description A structure to store isoform information.
struct Isoform
{
  string mature_id;
  string mature_seq;
  string seed_seq;
  string template_seq;
  string read_seq;
  int tpm;
  int indel_5p;
  int indel_3p;
  int snp_5p;
  int snp_3p;
  int snp_seed;

  Isoform(const string &mature_id, const string &mature_seq,
          const string &seed_seq, const string &template_seq,
          const string &read_seq, int tpm, int indel_5p,
          int indel_3p, int snp_5p, int snp_3p, int snp_seed)
      : mature_id(mature_id), mature_seq(mature_seq), seed_seq(seed_seq),
        template_seq(template_seq), read_seq(read_seq), tpm(tpm),
        indel_5p(indel_5p), indel_3p(indel_3p), snp_5p(snp_5p), snp_3p(snp_3p),
        snp_seed(snp_seed) {}

  int get_dist()
  {
    int dist = indel_5p + indel_3p + snp_5p + snp_3p + snp_seed;
    return dist;
  }
};

// @title Mirna Structure
// @description A structure to store miRNA information.
struct Mirna
{
  string mature_id;
  string mature_seq;
  string seed_seq;
  string template_seq;
  string flank_5p_seq;
  string flank_3p_seq;

  Mirna(const string &mature_id, const string &mature_seq, const string &seed_seq,
        const string &template_seq, const string &flank_5p_seq,
        const string &flank_3p_seq)
      : mature_id(mature_id),
        mature_seq(mature_seq),
        seed_seq(seed_seq),
        template_seq(template_seq),
        flank_5p_seq(flank_5p_seq),
        flank_3p_seq(flank_3p_seq) {}
};

// @title Calculate Edit Distance
// @description Calculate the edit distance between two strings.
// @param s1 First string
// @param s2 Second string
// @return The edit distance between s1 and s2
int calc_edit_dist(const string &s1, const string &s2)
{
  int m = s1.size();
  int n = s2.size();

  vector<vector<int>> dp(m + 1, vector<int>(n + 1));

  for (int i = 0; i <= m; i++)
    dp[i][0] = i;
  for (int j = 0; j <= n; j++)
    dp[0][j] = j;

  for (int i = 1; i <= m; i++)
  {
    for (int j = 1; j <= n; j++)
    {
      if (s1[i - 1] == s2[j - 1])
        dp[i][j] = dp[i - 1][j - 1];
      else
        dp[i][j] = 1 + min({dp[i][j - 1],
                            dp[i - 1][j],
                            dp[i - 1][j - 1]});
    }
  }

  return dp[m][n];
}

pair<int, int> find_seed(const string &read_seq, const string &seed_seq,
                         int max_snp)
{
  if (read_seq.length() < 16)
    stop("The length of read: %s is too short", read_seq);

  if (seed_seq.length() > 16 || seed_seq.length() < 8)
    stop("The length of seed: %s is shoule be between 8 ant 16 nt", seed_seq);

  uint64_t mask = 0xFFFFFFFFFFFFFFFF;

  if (seed_seq.length() < 16)
    mask >>= (16 - seed_seq.length()) * 4;

  vector<uint64_t> packed_reads = pack_read(read_seq, seed_seq.length(), mask);
  uint64_t packed_seed = pack_seed(seed_seq);

  int max_diff = max_snp * 2;
  int min_diff = 64;
  int index = -1;

  for (int i = 0; i < packed_reads.size(); i++)
  {
    int diff = bit_count(packed_reads[i] ^ packed_seed);
    if (diff <= max_diff && diff < min_diff)
    {
      index = i;
      min_diff = diff;
    }
    if (diff == 0)
      break;
  }

  pair<int, int> ret;
  ret.first = index;
  ret.second = min_diff / 2;

  return ret;
}

pair<int, int> calc_dist_5p(const string &read_5p, const string &mature_5p,
                            const string &flank_5p)
{
  int offset = read_5p.length() - mature_5p.length();
  int indel = abs(offset);
  int snp = calc_snp_5p(read_5p, mature_5p);

  if (offset > 0)
  {
    string read_5p_1 = read_5p.substr(0, offset);
    string read_5p_2 = read_5p.substr(offset);

    if (read_5p_1.length() > flank_5p.length())
    {
      indel = -1;
      snp = -1;
    }
    else
    {
      snp += calc_snp_5p(read_5p_1, flank_5p);
    }
  }

  pair<int, int> ret;
  ret.first = indel;
  ret.second = snp;

  return ret;
}

pair<int, int> calc_dist_3p(const string &read_3p, const string &mature_3p)
{
  int offset = read_3p.length() - mature_3p.length();
  int indel = abs(offset);
  int snp = calc_snp_3p(read_3p, mature_3p);
  pair<int, int> ret;
  ret.first = indel;
  ret.second = snp;

  return ret;
}

vector<Isoform *> detect_one_seq(const string &read_seq, double tpm,
                                 const vector<Mirna *> &mirnas,
                                 int max_indel_5p, int max_indel_3p,
                                 int max_snp_5p, int max_snp_3p,
                                 int max_snp_seed, int max_snp)
{
  vector<Isoform *> hits;
  vector<Isoform *> reduced_hits;

  for (Mirna *mirna : mirnas)
  {
    string mature_id = mirna->mature_id;
    string mature_seq = mirna->mature_seq;
    string seed_seq = mirna->seed_seq;
    string template_seq = mirna->template_seq;
    string flank_5p_seq = mirna->flank_5p_seq;
    string flank_3p_seq = mirna->flank_3p_seq;

    auto seed = find_seed(read_seq, seed_seq, max_snp_seed);
    int read_seq_index_5p = seed.first;
    int snp_seed = seed.second;
    if (read_seq_index_5p == -1)
      continue;

    size_t mature_index_5p = mature_seq.find(seed_seq);
    if (mature_index_5p == string::npos)
      stop("Seed: %s is not found in reference", seed_seq);

    size_t mature_index_3p = mature_index_5p + seed_seq.length();
    size_t read_seq_index_3p = read_seq_index_5p + seed_seq.length();
    string mature_5p = mature_seq.substr(0, mature_index_5p);
    string mature_3p = mature_seq.substr(mature_index_3p);
    string read_seq_5p = read_seq.substr(0, read_seq_index_5p);
    string read_seq_3p = read_seq.substr(read_seq_index_3p);

    auto dist_5p = calc_dist_5p(read_seq_5p, mature_5p, flank_5p_seq);
    int indel_5p = dist_5p.first;
    int snp_5p = dist_5p.second;
    if (indel_5p == -1)
      continue;

    auto dist_3p = calc_dist_3p(read_seq_3p, mature_3p);
    int indel_3p = dist_3p.first;
    int snp_3p = dist_3p.second;

    int snp = snp_5p + snp_3p + snp_seed;

    if (indel_5p <= max_indel_5p && indel_3p <= max_indel_3p &&
        snp_5p <= max_snp_5p && snp_3p <= max_snp_3p && snp_seed <= max_snp_seed
        && snp <= max_snp)
    {
      hits.push_back(new Isoform(mature_id, mature_seq, seed_seq, template_seq,
                                 read_seq, tpm, indel_5p, indel_3p,
                                 snp_5p, snp_3p, snp_seed));
    }
  }

  if (!hits.empty())
  {
    int min_dist = hits[0]->get_dist();

    for (int i = 1; i < hits.size(); i++)
    {
      Isoform *hit = hits[i];
      if (hit->get_dist() < min_dist)
        min_dist = hit->get_dist();
    }

    for (int i = 0; i < hits.size(); i++)
    {
      if (hits[i]->get_dist() == min_dist)
        reduced_hits.push_back(hits[i]);
    }
  }

  return reduced_hits;
}

vector<Mirna *> load_mirs(const vector<string> &mature_ids,
                          const vector<string> &mature_seqs,
                          const vector<string> &seed_seqs,
                          const vector<string> &template_seqs,
                          const vector<string> &flank_5p_seqs,
                          const vector<string> &flank_3p_seqs)
{
  vector<Mirna *> mirnas;
  for (int i = 0; i < mature_ids.size(); i++)
  {
    Mirna *mirna = new Mirna(mature_ids[i], mature_seqs[i], seed_seqs[i],
                             template_seqs[i], flank_5p_seqs[i], flank_3p_seqs[i]);
    mirnas.push_back(mirna);
  }
  return mirnas;
}

//' @title Find Isoforms
//' @description Find isoforms for given miRNAs and reads.
//' @param mirnas A DataFrame containing miRNA information. The DataFrame should
//'   include the following columns:
//'   \itemize{
//'     \item mature_id: Unique identifier of the mature miRNA.
//'     \item mature_seq: Sequence of the mature miRNA.
//'     \item seed_seq: Seed sequence of the mature miRNA.
//'     \item template_seq: Template sequence including flanking regions.
//'     \item flank_5p_seq: 5' flanking sequence of the mature miRNA.
//'     \item flank_3p_seq: 3' flanking sequence of the mature miRNA.
//'   }
//' @param reads A DataFrame containing read information. The DataFrame should
//'   include the following columns:
//'   \itemize{
//'     \item read_seq: Sequence of the read.
//'     \item tpm: TPM of the read.
//'   }
//' @param max_indel_5p An integer specifying the maximum allowed indels at the
//'   5' end. Must be non-negative.
//' @param max_indel_3p An integer specifying the maximum allowed indels at the
//'   3' end. Must be non-negative.
//' @param max_snp_5p An integer specifying the maximum allowed mismatches at
//'   the 5' end. Must be non-negative.
//' @param max_snp_3p An integer specifying the maximum allowed mismatches at
//'   the 3' end. Must be non-negative.
//' @param max_snp_seed An integer specifying the maximum allowed mismatches in
//'   the seed region. Must be non-negative.
//' @param max_snp An integer specifying the maximum allowed mismatches. Must be
//'   non-negative.
//' @return A DataFrame containing the detected isoforms with the following columns:
//'   \itemize{
//'     \item mature_id: Reference miRNA ID.
//'     \item mature_seq: Reference miRNA sequence.
//'     \item seed_seq: Seed sequence of the mature miRNA.
//'     \item template_seq: Template sequence including flanking regions.
//'     \item read_seq: Read sequence matched to the reference.
//'     \item read_num: Abundance of the read.
//'     \item tpm: TPM of the read.
//'     \item indel_5p: Indels at the 5' end.
//'     \item indel_3p: Indels at the 3' end.
//'     \item snp_5p: Mismatches at the 5' end.
//'     \item snp_3p: Mismatches at the 3' end.
//'     \item snp_seed: Mismatches in the seed.
//'     \item dist: Total edit distance.
//'   }
//' @details
//' This function compares read sequences to reference miRNAs and identifies
//' isoforms based on the specified thresholds for indels and mismatches. It calculates the edit distance between the read sequences and the reference miRNAs, and returns a DataFrame containing the detected isoforms that meet the criteria.
//' @examples
//' \dontrun{
//' }
//' @export
// [[Rcpp::export]]
DataFrame find_isoforms(DataFrame mirnas, DataFrame reads,
                        int max_indel_5p, int max_indel_3p, int max_snp_5p,
                        int max_snp_3p, int max_snp_seed, int max_snp)
{
  if (max_indel_5p < 0 || max_indel_3p < 0 || max_snp_5p < 0 ||
      max_snp_3p < 0 || max_snp_seed < 0 || max_snp < 0)
  {
    stop("All indel and SNP thresholds must be non-negative.");
  }

  vector<string> mature_ids = as<vector<string>>(mirnas["mature_id"]);
  vector<string> mature_seqs = as<vector<string>>(mirnas["mature_seq"]);
  vector<string> seed_seqs = as<vector<string>>(mirnas["seed_seq"]);
  vector<string> template_seqs = as<vector<string>>(mirnas["template_seq"]);
  vector<string> flank_5p_seqs = as<vector<string>>(mirnas["flank_5p_seq"]);
  vector<string> flank_3p_seqs = as<vector<string>>(mirnas["flank_3p_seq"]);
  vector<string> read_seqs = as<vector<string>>(reads["read_seq"]);
  vector<double> tpms = as<vector<double>>(reads["tpm"]);

  vector<Mirna *> matures = load_mirs(mature_ids, mature_seqs, seed_seqs,
                                      template_seqs, flank_5p_seqs, flank_3p_seqs);

  vector<Isoform *> isoforms;

  for (size_t i = 0; i < read_seqs.size(); ++i)
  {
    const string &read_seq = read_seqs[i];
    double tpm = tpms[i];

    vector<Isoform *> hits = detect_one_seq(
        read_seq, tpm, matures, max_indel_5p, max_indel_3p, max_snp_5p,
        max_snp_3p, max_snp_seed, max_snp);

    isoforms.insert(isoforms.end(), hits.begin(), hits.end());
  }

  vector<string> result_mature_ids;
  vector<string> result_mature_seqs;
  vector<string> result_seed_seqs;
  vector<string> result_template_seqs;
  vector<string> result_read_seqs;
  vector<int> result_tpms;
  vector<int> result_indels_5p;
  vector<int> result_indels_3p;
  vector<int> result_snps_5p;
  vector<int> result_snps_3p;
  vector<int> result_snps_seed;
  vector<int> result_dists;

  for (Isoform *isoform : isoforms)
  {
    result_mature_ids.push_back(isoform->mature_id);
    result_mature_seqs.push_back(isoform->mature_seq);
    result_seed_seqs.push_back(isoform->seed_seq);
    result_template_seqs.push_back(isoform->template_seq);
    result_read_seqs.push_back(isoform->read_seq);
    result_tpms.push_back(isoform->tpm);
    result_indels_5p.push_back(isoform->indel_5p);
    result_indels_3p.push_back(isoform->indel_3p);
    result_snps_5p.push_back(isoform->snp_5p);
    result_snps_3p.push_back(isoform->snp_3p);
    result_snps_seed.push_back(isoform->snp_seed);
    result_dists.push_back(isoform->get_dist());

    delete isoform;
  }

  for (Mirna *mirna : matures)
  {
    delete mirna;
  }

  return DataFrame::create(
      _["mature_id"] = result_mature_ids,
      _["mature_seq"] = result_mature_seqs,
      _["seed_seq"] = result_seed_seqs,
      _["template_seq"] = result_template_seqs,
      _["read_seq"] = result_read_seqs,
      _["tpm"] = result_tpms,
      _["indel_5p"] = result_indels_5p,
      _["indel_3p"] = result_indels_3p,
      _["snp_5p"] = result_snps_5p,
      _["snp_3p"] = result_snps_3p,
      _["snp_seed"] = result_snps_seed,
      _["dist"] = result_dists);
}
