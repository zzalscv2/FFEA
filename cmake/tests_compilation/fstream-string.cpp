#include <string>
#include <fstream>

int main () {

  std::ifstream fin;
  std::string k = "whatever";
  fin.open(k);
  fin.close();
  return 0; 

}
