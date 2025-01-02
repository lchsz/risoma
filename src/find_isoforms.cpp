#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


struct Isoform
{
  string mature_id;
  string mature_seq;
  string template_seq;
  string read_seq;
  int read_num;
  int dist;

  Isoform(const string &mature_id, const string &mature_seq,
          const string &template_seq, const string &read_seq, int read_num,
          int dist)
      : mature_id(mature_id), mature_seq(mature_seq), template_seq(template_seq),
        read_seq(read_seq), read_num(read_num), dist(dist) {}
};


struct Mirna
{
  string mature_id;
  string mature_seq;
  string seed_seq;
  string template_seq;

  Mirna(const string &mature_id, const string &mature_seq,
        const string &seed_seq, const string &template_seq)
      : mature_id(mature_id), mature_seq(mature_seq), seed_seq(seed_seq),
        template_seq(template_seq){}
};


//' Find the minimum number of edits to convert ‘s1‘ into ‘s2‘.
//'
//' @param s1 first string
//' @param s2 second string
//' @return edit distance
// [[Rcpp::export]]
int calc_edit_dist(const std::string &s1, const std::string &s2)
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


vector<Isoform *> detect_one_seq(const string &read_seq, int read_num,
                               const vector<Mirna *> &mirnas,
                               int max_ed_5p, int max_ed_3p, const string &type)
{
  vector<Isoform *> hits;
  vector<Isoform *> reduced_hits;

  for (Mirna *mirna : mirnas)
  {
    string mature_id = mirna->mature_id;
    string mature_seq = mirna->mature_seq;
    string template_seq = mirna->template_seq;
    string seed_seq = mirna->seed_seq;

    size_t read_seq_index_5p = read_seq.find(seed_seq);
    if (read_seq_index_5p == string::npos)
      continue;

    size_t mature_index_5p = mature_seq.find(seed_seq);
    if (mature_index_5p == string::npos)
    {
      stop("Seed: %s is not found in consensus", seed_seq);
    }

    size_t mature_index_3p = mature_index_5p + seed_seq.size();
    size_t read_seq_index_3p = read_seq_index_5p + seed_seq.size();
    string mature_5p = mature_seq.substr(0, mature_index_5p);
    string mature_3p = mature_seq.substr(mature_index_3p);
    string read_seq_5p = read_seq.substr(0, read_seq_index_5p);
    string read_seq_3p = read_seq.substr(read_seq_index_3p);

    int dist_5p = calc_edit_dist(read_seq_5p, mature_5p);
    int dist_3p = calc_edit_dist(read_seq_3p, mature_3p);

    bool is_isoform = false;
    if (type == "all") {
      if (dist_5p <= max_ed_5p && dist_3p <= max_ed_3p) {
        is_isoform = true;
      }
    } else if (type == "p5") {
      if (dist_5p < max_ed_5p && dist_3p == 0) {
        is_isoform = true;
      }
    } else if (type == "p3") {
      if (dist_3p <= max_ed_3p && dist_5p == 0) {
        is_isoform = true;
      }
    }

    if (is_isoform) {
      hits.push_back(new Isoform(mature_id, mature_seq, template_seq, read_seq,
                                 read_num, dist_5p + dist_3p));
    }
  }

  if (!hits.empty())
  {
    int min_dist = hits[0]->dist;

    for (int i = 1; i < hits.size(); i++)
    {
      Isoform *hit = hits[i];
      if (hit->dist < min_dist)
        min_dist = hit->dist;
    }

    for (int i = 0; i < hits.size(); i++)
    {
      if (hits[i]->dist == min_dist)
        reduced_hits.push_back(hits[i]);
    }
  }

  return reduced_hits;
}


vector<Mirna *> load_mirs(const vector<string> &mature_ids,
                         const vector<string> &mature_seqs,
                         const vector<string> &seed_seqs,
                         const vector<string> &template_seqs)
{
  vector<Mirna *> mirnas;
  for (int i = 0; i < mature_ids.size(); i++)
  {
    Mirna *mirna = new Mirna(mature_ids[i], mature_seqs[i], seed_seqs[i],
                             template_seqs[i]);
    mirnas.push_back(mirna);
  }
  return mirnas;
}


//' Find isoform from reads based on reference miRNAs
//'
//' @param mirnas a data.frame storing mature miRNA ID (matrue_ID),
//'   mature miRNA sequence (mature_seq), mature start (mature_start),
//'   seed on miRNA (seed_seq), precursor ID (pre_ID) and precursor sequence
//'   (pre_seq)
//' @param reads a data.frame containing sequence and it's amount
//' @param maxEd5p the maximum distance between 5’ region of reads with
//'   reference miRNA
//' @param maxEd3p the maximum distance between 3' region of reads with
//'   reference miRNA
//' @return a data.frame with columns: mature_ID, read_seq, read_num and dist
// [[Rcpp::export]]
DataFrame find_isoforms(DataFrame mirnas, DataFrame reads,
                       int max_ed_5p, int max_ed_3p, std::string type)
{
  vector<string> mature_ids = as<vector<string>>(
      as<CharacterVector>(mirnas["mature_id"]));
  vector<string> mature_seqs = as<vector<string>>(
      as<CharacterVector>(mirnas["mature_seq"]));
  vector<string> seed_seqs = as<vector<string>>(
      as<CharacterVector>(mirnas["seed_seq"]));
  vector<string> template_seqs = as<vector<string>>(
      as<CharacterVector>(mirnas["template_seq"]));

  vector<string> read_seqs = as<vector<string>>(
      as<CharacterVector>(reads["read_seq"]));
  vector<int> read_nums = as<vector<int>>(
      as<IntegerVector>(reads["read_num"]));

  vector<Mirna *> matures = load_mirs(mature_ids, mature_seqs, seed_seqs,
                                      template_seqs);

  vector<Isoform *> isoforms;

  for (int i = 0; i < read_seqs.size(); i++)
  {
    string read_seq = read_seqs[i];
    int read_num = read_nums[i];
    vector<Isoform *> hits = detect_one_seq(read_seq, read_num, matures,
                                          max_ed_5p, max_ed_3p, type);
    isoforms.insert(isoforms.end(), hits.begin(), hits.end());
  }

  vector<string> result_mature_ids;
  vector<string> result_mature_seqs;
  vector<string> result_template_seqs;
  vector<string> result_read_seqs;
  vector<int> result_read_nums;
  vector<int> result_dists;

  for (Isoform *isoform : isoforms)
  {
    result_mature_ids.push_back(isoform->mature_id);
    result_mature_seqs.push_back(isoform->mature_seq);
    result_template_seqs.push_back(isoform->template_seq);
    result_read_seqs.push_back(isoform->read_seq);
    result_read_nums.push_back(isoform->read_num);
    result_dists.push_back(isoform->dist);
  }

  DataFrame result = DataFrame::create(
      _["mature_id"] = result_mature_ids,
      _["mature_seq"] = result_mature_seqs,
      _["template_seq"] = result_template_seqs,
      _["read_seq"] = result_read_seqs,
      _["read_num"] = result_read_nums,
      _["dist"] = result_dists);

  return result;
}
