#include <iostream>

typedef double scalar;
typedef scalar arr3[3];


class vector3 {
public:
    arr3 data;
    scalar& x = data[0];
    scalar& y = data[1];
    scalar& z = data[2];
    scalar& operator [](std::size_t i) { return data[i]; }
    void assign( scalar t0, scalar t1, scalar t2 )
    { data[0] = t0; data[1] = t1; data[2] = t2; }
};

int main() {

  vector3 v;
  v.assign(0, 1, 2);
  std::cout << v[0] << ", " << v[1] << ", " << v[2] << std::endl; 

  return 0;

}
