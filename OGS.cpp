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
void writeMatrixToFile(const char *, mat_RR &);

mat_ZZ B;
mat_RR Bp, mi;
int n;

int main(int argc, char const *argv[]) {
  readMatrixFromFile(argv[1], B);
  performGSO(Bp, mi);
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

void writeMatrixToFile(const char *str, mat_RR &M) {
  std::ofstream writer(str);
  writer << n << "\n";
  for(int i = 0; i < n; i++) writer << M[i] << "\n";
  writer.close();
}
