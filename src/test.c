#include <fstream.h>

void main () {

  ifstream f1;
  ofstream f2;
  f1.open("scores.96");
  f2.open("final.96");

  int s1, s2, s3;
  float w1, w2, w3;

  f1 >> s1 >> w1;
  f1 >> s2 >> w2;
  f1 >> s3 >> w3;
  f2 << (s1*w1+s2*w2+s3*w3);

}
