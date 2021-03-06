#include "lattice/rpbkz.hpp"

void ExtractLWEChallengeFile(mat_ZZ& A, vec_ZZ& b, double& sigma2, int& q,
                             std::string filename) {
  ifstream ff;
  ff.open(filename.c_str(), ios_base::in);
  int n, m;
  ff >> q;
  cout << "q=" << q << endl;
  ff >> sigma2;
  cout << "sigma2=" << sigma2 << endl;
  ff >> b;
  cout << "dim(b)=" << b.length() << endl;
  ff >> A;
  cout << "size(A)=" << A.NumRows() << "," << A.NumCols() << endl;
}

bool solve_lwe(char* fn, double coe, int core) {
  std::string lwefile, outfile, basfile;
  lwefile = "data/" + string(fn) + ".txt";
  outfile = "sol/" + string(fn) + ".txt";
  if (FileExists(lwefile) == false) {
    cout << "file(" << lwefile << ") does not exist" << endl;
    cout << "skip this test" << endl;
    return true;
  }

  mat_ZZ lweA;
  vec_ZZ lweb;
  double lwesigma2;
  int lweq;

  ExtractLWEChallengeFile(lweA, lweb, lwesigma2, lweq, lwefile);

  int n = lweA.NumCols();
  int m = lweA.NumRows();
  m = std::min(m, int(n * 2));
  ZZ M = to_ZZ(1);
  cout << "subdim=" << m << endl;
  double sigma = sqrt(lwesigma2);
  cout << "sigma=" << sigma << endl;
  bkzfloat maxnorm;
  maxnorm = coe * sqrt(m) * sigma;
  cout << "maxnorm of e=" << maxnorm << endl;

  vec_ZZ target;
  target.SetLength(m);
  for (int i = 0; i < m; i++) target[i] = lweb[i];

  LatticeBasis<long double> A, B;
  A.resize(m + n, m);
  for (int i = 0; i < m; i++) A.L[i][i] = lweq;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) A.L[i + m][j] = lweA[j][i];
  }
  cout << "Find basis vectors for A" << endl;
  ::BigLLL(A.L, 0, 0.999, VL1);
  cout << "Apply PBKZ" << endl;
  B.resize(m, m);
  for (int i = 0; i < m; i++) B.L[i] = A.L[i];
  ProgressiveBKZ(B, 0, 45, VL1, "ignoreflat");
  ofstream bff;
  bff.open(basfile.c_str(), ofstream::out);
  bff << B.L << endl;
  cout << "Find close vector" << endl;
  mat_ZZ VV;
  std::ostringstream stringStream;
  stringStream << "parallel=" << core;
  string otheroption = stringStream.str();
  VV = ENUMCV(B, target, maxnorm, 0.1, enum_mode_all_vectors, 0, VL3,
              otheroption);
  vec_ZZ candidate;
  int can = -1;  // candidate index
  for (int i = 0; i < VV.NumRows(); i++) {
    candidate = VV[i] - target;
    cout << "error candidate[" << i + 1 << "]: " << candidate << endl;
    cout << "norm[" << i + 1 << "]=" << LengthOf(candidate) << endl;
    if (LengthOf(candidate) < maxnorm) {
      can = i;
      break;
    }
  }
  if (can == -1) return false;

  // From error to secret by gaussian elimination
  cout << "Extract secret vector from error candidate by gaussian elimitation" << endl;
  for (int i = 0; i < n; i++) {
    // normalize i-th row
    ZZ a;
    int shift = 0;
    while (GCD(lweA[i + shift][i], to_ZZ(lweq)) != 1 && i + shift < m) shift++;
    if (i + shift == m) cout << "ERROR!" << endl;
    swap(lweA[i], lweA[i + shift]);
    swap(VV[can][i], VV[can][i + shift]);
    a = lweA[i][i];
    ZZ ai = InvMod(a, to_ZZ(lweq));
    for (int j = 0; j < n; j++) {
      lweA[i][j] *= ai;
      lweA[i][j] %= to_ZZ(lweq);
      if (lweA[i][j] < 0) lweA[i][j] += lweq;
    }
    VV[can][i] *= ai;
    VV[can][i] %= to_ZZ(lweq);
    if (VV[can][i] < 0) VV[can][i] += lweq;

    for (int j = 0; j < m; j++) {
      if (j != i) {
        VV[can][j] -= lweA[j][i] * VV[can][i];
        lweA[j] -= lweA[j][i] * lweA[i];
      }
      for (int k = 0; k < n; k++) {
        lweA[j][k] %= to_ZZ(lweq);
        if (lweA[j][k] < 0) lweA[j][k] += lweq;
      }
      VV[can][j] %= to_ZZ(lweq);
      if (VV[can][j] < 0) VV[can][j] += lweq;
    }
  }

  ofstream fff;
  fff.open(outfile.c_str(), ofstream::out | ofstream::app);

  cout << "final solution s=[";
  for (int i = 0; i < n; i++) cout << VV[can][i] << " ";
  cout << "]" << endl;

  fff << "[";
  for (int i = 0; i < n; i++) fff << VV[can][i] << " ";
  fff << "]" << endl;

  cout << "VV=[";
  for (int i = 0; i < m; i++) cout << VV[can][i] << " ";
  cout << "]" << endl;

  return true;
}

int main(int argc, char** argv) {
  solve_lwe(argv[1], atof(argv[2]), atoi(argv[3]));
  exit(0);
}
