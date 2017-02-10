#include <future>

unsigned int test(int i)
{

  return i;

}

int main() {

  std::future<unsigned int> f1;
  f1 = std::async(std::launch::async,test,3);
  f1.wait();
  unsigned int j = f1.get();
  return (j == 3) ? 0: 1;

}
  

