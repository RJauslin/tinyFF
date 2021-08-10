#include <RcppThread.h>
// [[Rcpp::export]]
std::vector<double> test2()
{
  RcppThread::ThreadPool pool;
  auto task = [] (int i) {
    double result(0.0);
    result = result + i;
    return result;
  };
  std::vector<std::future<double>> futures(10);
  std::vector<double> results(10);
  for (unsigned int i = 0; i < 10; ++i)
    futures[i] = pool.pushReturn(task, i);
  for (unsigned int i = 0; i < 10; ++i)
    results[i] = futures[i].get();
  pool.join();
  return(results);
}