#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/RR.h>

using namespace NTL;

void readMatrixFromFile(const char *, mat_ZZ &);
void computeHB(RR &, mat_ZZ &);

mat_ZZ B;
int n;
RR HB;

int main(int argc, char const *argv[]) {
  readMatrixFromFile(argv[1], B);
  computeHB(HB, B);
  std::cout << HB << "\n";

  return 0;
}

void readMatrixFromFile(const char *str, mat_ZZ &B) {
  std::ifstream reader(str);
  reader >> n;
  B.SetDims(n, n);
  for(int i = 0; i < n; i++) reader >> B[i];
  reader.close();
}

void computeHB(RR &HB, mat_ZZ &B) {
  ZZ num;
  RR den;
  vec_ZZ inners;
  inners.SetLength(n);

  determinant(num, B, 0);
  abs(num, num);
  den = 1;
  for(int i = 0; i < n; i++) {
    InnerProduct(inners[i], B[i], B[i]);
    den *= sqrt(conv<RR>(inners[i]));
  }
  div(HB, conv<RR>(num), den);
}
