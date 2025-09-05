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
int bw_x = 4;

uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

void test_diff(uint64_t *x) {
  uint64_t *y = new uint64_t[dim];

  uint64_t comm_start = iopack->get_comm();
  INIT_TIMER;
  START_TIMER;
  math->diff(dim, x, y, bw_x);
  STOP_TIMER("Diff");
  double total_comm = (iopack->get_comm() - comm_start) / (1.0 * (1ULL << 20)); // In MB
  cout << "Total Communication for Diff: " << total_comm << " mb" << endl;

  if (party == ALICE) {
    iopack->io->send_data(x, dim * sizeof(uint64_t));
    iopack->io->send_data(y, dim * sizeof(uint64_t));
  } else { // party == BOB
    uint64_t *x0 = new uint64_t[dim];
    uint64_t *y0 = new uint64_t[dim];
    iopack->io->recv_data(x0, dim * sizeof(uint64_t));
    iopack->io->recv_data(y0, dim * sizeof(uint64_t));

    for (int i = 0; i < dim; i++) {
      uint64_t D = (x[i] - x0[i]) & mask_x;
      uint64_t Y = (y[i] + y0[i]) & mask_x;
      assert(Y == D);
    }
    cout << "Diff Tests Passed" << endl;

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

  test_diff(x);

  delete[] x;
  delete math;
}
