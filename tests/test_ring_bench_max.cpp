#include "Math/math-functions.h"
#include "utils/emp-tool.h"
#include <iostream>

using namespace sci;
using namespace std;

int party, port = 32000;
string address = "127.0.0.1";
IOPack *iopack;
OTPack *otpack;
MathFunctions *math;

int dim = 1ULL << 16;
int bw_x = 32;

uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

uint64_t find_max(uint64_t *x, int dim) {
  if (dim == 0) {
    throw std::invalid_argument("find_max: dim must be greater than 0");
  }
  uint64_t max = x[0];
  for (int i = 1; i < dim; i++) {
    if (x[i] > max) {
      max = x[i];
    }
  }
  return max;
}

void test_max(uint64_t *x) {
  uint64_t *y = new uint64_t[1];

  uint64_t comm_start = math->iopack->get_comm();
  INIT_TIMER;
  START_TIMER;
  math->max(dim, x, y, bw_x);
  // math->max_naive(dim, x, y, bw_x);
  STOP_TIMER("Max:");
  uint64_t comm_end = math->iopack->get_comm();
  double comm = (comm_end - comm_start) / (1.0 * (1ULL << 20));
  cout << "Total Communication for Max: " << comm << " mb" << endl;

  if (party == ALICE) {
    iopack->io->send_data(x, dim * sizeof(uint64_t));
    iopack->io->send_data(y, sizeof(uint64_t));
  } else { // party == BOB
    uint64_t *x0 = new uint64_t[dim];
    uint64_t *y0 = new uint64_t[1];
    iopack->io->recv_data(x0, dim * sizeof(uint64_t));
    iopack->io->recv_data(y0, sizeof(uint64_t));

    // Compute the absolute value of the differences
    for (int i = 0; i < dim; i++) {
      x[i] = (x[i] + x0[i]) & mask_x;
      int64_t X = signed_val(x[i], bw_x);
      if (X < 0)
        X = -X;
      x[i] = unsigned_val(X, bw_x);
    }

    uint64_t max = find_max(x, dim);
    
    uint64_t Y = (y[0] + y0[0]) & mask_x;
    
    assert(Y == max);

    cout << "Max Tests Passed" << endl;

    delete[] x0;
    delete[] y0;
  }
  delete[] y;
}

int main(int argc, char **argv) {
  ArgMapping amap;
  amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
  amap.arg("p", port, "Port Number");
  amap.arg("ip", address, "IP Address of server (ALICE)");

  amap.parse(argc, argv);

  iopack = new IOPack(party, port, address);
  otpack = new OTPack(iopack, party);

  math = new MathFunctions(party, iopack, otpack);

  PRG128 prg;

  uint64_t *x = new uint64_t[dim];

  prg.random_data(x, dim * sizeof(uint64_t));

  for (int i = 0; i < dim; i++) {
    x[i] &= mask_x;
  }
  math->ABS(dim, x, x, bw_x);

  test_max(x);

  delete[] x;
  delete math;
}
