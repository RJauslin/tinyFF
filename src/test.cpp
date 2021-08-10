#include <RcppThread.h>
// [[Rcpp::export]]
void greet()
{
  auto job = [] () {
    for (size_t i = 0; i < 100; ++i)
      RcppThread::Rcout << "Hi!" << std::endl;
  };
  std::thread t1(job);
  std::thread t2(job);
  t1.join();
  t2.join();
  RcppThread::Rcout << "";
}
