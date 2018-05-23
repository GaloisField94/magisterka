#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/RR.h>
#include <NTL/vec_RR.h>
#include <NTL/mat_RR.h>

using namespace NTL;

void readMatrixFromFile(const char *, mat_ZZ &);
void performGSO(mat_RR &, mat_RR &);
void performSR(mat_ZZ &, mat_RR &);
bool LovashCond(mat_RR &, mat_RR &, RR, int);
void performBasicLLL(mat_ZZ &, mat_RR &, mat_RR &, RR);
void writeMatrixToFile(const char *, mat_RR &);
void writeMatrixToFile(const char *, mat_ZZ &);

mat_ZZ B;
mat_RR Bp, mi;
RR delta = RR(0.75);
int n;

int main(int argc, char const *argv[]) {
  readMatrixFromFile(argv[1], B);
  performBasicLLL(B, Bp, mi, delta);
  writeMatrixToFile("Bblll.txt", B);
  writeMatrixToFile("Bp.txt", Bp);
  writeMatrixToFile("mi.txt", mi);

  return 0;
}

void readMatrixFromFile(const char *str, mat_ZZ &B) {
  std::ifstream reader(str);
  reader >> n;
  B.SetDims(n, n);
  for(int i = 0; i < n; i++) reader >> B[i];
  reader.close();
}

void performGSO(mat_RR &Bp, mat_RR &mi) {
  RR mi_nom, mi_den;
  vec_RR sum;
  sum.SetLength(n);

  Bp.SetDims(n, n);
  Bp[0] = conv<vec_RR>(B[0]);
  mi.SetDims(n, n);
  for(int i = 0; i < n; i++) mi[i][i] = 1;

  for(int i = 1; i < n; i++) {
    for(int j = 0; j < i; j++) {
      InnerProduct(mi_nom, conv<vec_RR>(B[i]), Bp[j]);
      InnerProduct(mi_den, Bp[j], Bp[j]);
      mi[i][j] = mi_nom / mi_den;
      sum += mi[i][j] * Bp[j];
    }
    Bp[i] = conv<vec_RR>(B[i]) - sum;
    clear(sum);
  }
}

void performSR(mat_ZZ &B, mat_RR &mi) {
  ZZ tmp;

  for(int i = 1; i < n; i++)
    for(int j = i - 1; j >= 0; j--) {
      RoundToZZ(tmp, mi[i][j]);
      B[i] -= tmp * B[j];
      for(int k = 0; k <= j; k++) mi[i][k] -= conv<RR>(tmp) * mi[j][k];
    }
}

bool LovashCond(mat_RR &Bp, mat_RR &mi, RR delta, int i) {
  RR bi, bi1;
  bool satisfies;

  InnerProduct(bi, Bp[i], Bp[i]);
  InnerProduct(bi1, Bp[i - 1], Bp[i - 1]);
  std::cout << "||bi|| = " << bi << "\n";
  std::cout << "||bi1|| = " << bi1 << "\n";

  if(bi >= (delta - sqr(mi[i][i - 1])) * bi1) satisfies = true;
  else satisfies = false;

  return satisfies;
}

void performBasicLLL(mat_ZZ &B, mat_RR &Bp, mat_RR &mi, RR delta) {
  bool done;
  vec_ZZ tmp;
  tmp.SetLength(n);

  do {
    done = true;
    performGSO(Bp, mi);
    performSR(B, mi);
    std::cout << "Bp = " << Bp << "\n";
    std::cout << "mi = " << mi << "\n";
    for(int i = 1; i < n; i++) if(!LovashCond(Bp, mi, delta, i)) {
      done = false;
      tmp = B[i];
      B[i] = B[i - 1];
      B[i - 1] = tmp;
      break;
    }
  } while(!done);
}

void writeMatrixToFile(const char *str, mat_RR &M) {
  std::ofstream writer(str);
  writer << n << "\n";
  for(int i = 0; i < n; i++) writer << M[i] << "\n";
  writer.close();
}

void writeMatrixToFile(const char *str, mat_ZZ &M) {
  std::ofstream writer(str);
  writer << n << "\n";
  for(int i = 0; i < n; i++) writer << M[i] << "\n";
  writer.close();
}
