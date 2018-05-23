#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/RR.h>

using namespace NTL;

void readMatrixFromFile(const char *, mat_ZZ &);
void performLG(mat_ZZ &, mat_ZZ &);
void writeMatrixToFile(const char *, mat_ZZ &);

mat_ZZ B, B2;
int n;

int main(int argc, char const *argv[]) {
  readMatrixFromFile(argv[1], B);
  performLG(B2, B);
  writeMatrixToFile("LagGauB.txt", B2);

  return 0;
}

void readMatrixFromFile(const char *str, mat_ZZ &B) {
  std::ifstream reader(str);
  reader >> n;
  B.SetDims(n, n);
  for(int i = 0; i < n; i++) reader >> B[i];
  reader.close();
}

void performLG(mat_ZZ &B2, mat_ZZ &B) {
  ZZ num, den, q, b1_norm, b2_norm;
  vec_ZZ tmp;
  tmp.SetLength(n);
  B2.SetDims(n, n);

  InnerProduct(b1_norm, B[0], B[0]);
  InnerProduct(b2_norm, B[1], B[1]);
  if(b1_norm < b2_norm) {
    B2[0] = B[1];
    B2[1] = B[0];
  }
  else {
    B2[0] = B[0];
    B2[1] = B[1];
  }

  do {
    InnerProduct(num, B2[0], B2[1]);
    InnerProduct(den, B2[1], B2[1]);
    RoundToZZ(q, conv<RR>(num) / conv<RR>(den));
    tmp = B2[0] - q * B2[1];
    B2[0] = B2[1];
    B2[1] = tmp;
    InnerProduct(b1_norm, B2[0], B2[0]);
    InnerProduct(b2_norm, B2[1], B2[1]);
  } while(b1_norm > b2_norm);
}

void writeMatrixToFile(const char *str, mat_ZZ &M) {
  std::ofstream writer(str);
  writer << n << "\n";
  for(int i = 0; i < n; i++) writer << M[i] << "\n";
  writer.close();
}
