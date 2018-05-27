#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>

using namespace NTL;

void readVectorFromFile(const char *, vec_ZZ &);
void shiftRightVecZZ(vec_ZZ &);
void genBhq(mat_ZZ &);
void writeMatrixToFile(const char *, mat_ZZ &);

mat_ZZ Bhq;
vec_ZZ h;
uint16_t q;
int n;

int main(int argc, char const *argv[]) {
  readVectorFromFile(argv[1], h);
  genBhq(Bhq);
  writeMatrixToFile("Bhq.txt", Bhq);

  return 0;
}

void readVectorFromFile(const char *str, vec_ZZ &h) {
  std::ifstream reader(str);
  reader >> n;
  h.SetLength(2*n);
  reader >> q;
  reader >> h;
  reader.close();
}

void shiftRightVecZZ(vec_ZZ &v) {
  ZZ tmp;
  tmp = v[v.length() - 1];
  for(int i = v.length() - 1; i > 0; i--) v[i] = v[i-1];
  v[0] = tmp;
}

void genBhq(mat_ZZ &Bhq) {
  vec_ZZ zeroV, oneV, qV;

  zeroV.SetLength(n);
  clear(zeroV);

  oneV.SetLength(n);
  clear(oneV);
  oneV[0] = 1;

  qV.SetLength(n);
  clear(qV);
  qV[0] = q;

  Bhq.SetDims(2*n, 2*n);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      Bhq[i][j] = oneV[j];
      Bhq[i][j + n] = h[j];
    }
    shiftRightVecZZ(oneV);
    shiftRightVecZZ(h);
  }
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      Bhq[i + n][j] = zeroV[j];
      Bhq[i + n][j + n] = qV[j];
    }
    shiftRightVecZZ(qV);
  }
}

void writeMatrixToFile(const char *str, mat_ZZ &M) {
  std::ofstream writer(str);
  writer << 2*n << "\n";
  for(int i = 0; i < 2*n; i++) writer << M[i] << "\n";
  writer.close();
}
