#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

struct Isoform
{
  string matureId;
  string matureSeq;
  string readSeq;
  int readNum;
  int dist;

  Isoform(const string &matureId, const string &matureSeq, const string &readSeq,
          int readNum, int dist)
      : matureId(matureId), matureSeq(matureSeq), readSeq(readSeq),
        readNum(readNum), dist(dist) {}
};

struct Mirna
{
  string matureId;
  string matureSeq;
  string seedSeq;

  Mirna(const string &matureId, const string &matureSeq, const string &seedSeq)
      : matureId(matureId), matureSeq(matureSeq), seedSeq(seedSeq) {}
};

int findSeedInSeq(const string &seed, const string &seq)
{
  size_t index = seq.find(seed);
  if (index != string::npos)
  {
    return index;
  }
  else
  {
    return -1;
  }
}

//' Find the minimum number of edits to convert ‘s1‘ into ‘s2‘.
//'
//' @param s1 first string
//' @param s2 second string
//' @return edit distance
// [[Rcpp::export]]
int editDist(const std::string &s1, const std::string &s2)
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

int calcEditDist(const string &seq5p, const string &seq3p,
                 const string &consensus5p, const string &consensus3p,
                 int maxEd5p, int maxEd3p)
{
  int lev5p = editDist(seq5p, consensus5p);
  int lev3p = editDist(seq3p, consensus3p);
  int dist = lev3p + lev5p;
  if (maxEd5p != -1 && maxEd3p == -1)
  {
    if (lev5p > maxEd5p)
      dist = -5;
  }
  else if (maxEd5p == -1 && maxEd3p != -1)
  {
    if (lev3p > maxEd3p)
      dist = -3;
  }
  else if (maxEd5p != -1 && maxEd3p != -1)
  {
    if (lev5p > maxEd5p)
      dist = -5;
    if (lev3p > maxEd3p)
      dist = -3;
  }
  return dist;
}

vector<Isoform *> detectOneSeq(const string &readSeq, int readNum,
                               const vector<Mirna *> &mirnas,
                               int maxEd5p, int maxEd3p)
{
  vector<Isoform *> hits;
  vector<Isoform *> reducedHits;

  for (Mirna *mirna : mirnas)
  {
    string matureId = mirna->matureId;
    string matureSeq = mirna->matureSeq;
    string seedSeq = mirna->seedSeq;

    size_t readSeqIndex5p = readSeq.find(seedSeq);
    if (readSeqIndex5p == string::npos)
      continue;

    size_t matureIndex5p = matureSeq.find(seedSeq);
    if (matureIndex5p == string::npos)
    {
      stop("Seed: %s is not found in consensus", seedSeq);
    }

    size_t matureIndex3p = matureIndex5p + seedSeq.size();
    size_t readSeqIndex3p = readSeqIndex5p + seedSeq.size();
    string mature5p = matureSeq.substr(0, matureIndex5p);
    string mature3p = matureSeq.substr(matureIndex3p);
    string readSeq5p = readSeq.substr(0, readSeqIndex5p);
    string readSeq3p = readSeq.substr(readSeqIndex3p);
    int dist = calcEditDist(readSeq5p, readSeq3p, mature5p, mature3p,
                            maxEd5p, maxEd3p);
    if (dist == -3)
      continue;

    if (dist == -5)
      continue;

    hits.push_back(new Isoform(matureId, matureSeq, readSeq, readNum, dist));
  }

  if (!hits.empty())
  {
    int minDist = hits[0]->dist;

    for (int i = 1; i < hits.size(); i++)
    {
      Isoform *hit = hits[i];
      if (hit->dist < minDist)
        minDist = hit->dist;
    }

    for (int i = 0; i < hits.size(); i++)
    {
      if (hits[i]->dist == minDist)
        reducedHits.push_back(hits[i]);
    }
  }

  return reducedHits;
}

vector<Mirna *> loadMirs(const vector<string> &matureIds,
                         const vector<string> &matureSeqs,
                         const vector<string> &seedSeqs)
{
  vector<Mirna *> mirnas;
  for (int i = 0; i < matureIds.size(); i++)
  {
    Mirna *mirna = new Mirna(matureIds[i], matureSeqs[i], seedSeqs[i]);
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
DataFrame findIsoforms(DataFrame mirnas, DataFrame reads,
                       int maxEd5p, int maxEd3p)
{
  vector<string> matureIds = as<vector<string>>(
      as<CharacterVector>(mirnas["mature_ID"]));
  vector<string> matureSeqs = as<vector<string>>(
      as<CharacterVector>(mirnas["mature_seq"]));
  vector<string> seedSeqs = as<vector<string>>(
      as<CharacterVector>(mirnas["seed_seq"]));
  vector<string> preSeqs = as<vector<string>>(
      as<CharacterVector>(mirnas["pre_seq"]));

  vector<string> readSeqs = as<vector<string>>(
      as<CharacterVector>(reads["read_seq"]));
  vector<int> readNums = as<vector<int>>(
      as<IntegerVector>(reads["read_num"]));

  vector<Mirna *> matures = loadMirs(matureIds, matureSeqs, seedSeqs);

  vector<Isoform *> isoforms;

  for (int i = 0; i < readSeqs.size(); i++)
  {
    string readSeq = readSeqs[i];
    int readNum = readNums[i];
    vector<Isoform *> hits = detectOneSeq(readSeq, readNum, matures,
                                          maxEd5p, maxEd3p);
    isoforms.insert(isoforms.end(), hits.begin(), hits.end());
  }

  vector<string> resultMatureIds;
  vector<string> resultMatureSeqs;
  vector<string> resultReadSeqs;
  vector<int> resultReadNums;
  vector<int> resultDists;

  for (Isoform *isoform : isoforms)
  {
    resultMatureIds.push_back(isoform->matureId);
    resultMatureSeqs.push_back(isoform->matureSeq);
    resultReadSeqs.push_back(isoform->readSeq);
    resultReadNums.push_back(isoform->readNum);
    resultDists.push_back(isoform->dist);
  }

  DataFrame result = DataFrame::create(
      _["mature_ID"] = resultMatureIds,
      _["mature_seq"] = resultMatureSeqs,
      _["read_seq"] = resultReadSeqs,
      _["read_num"] = resultReadNums,
      _["dist"] = resultDists);
  
  return result;
}
