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

void test_adder(uint64_t *x, uint64_t *y) {
  uint64_t *z = new uint64_t[dim];

  uint64_t comm_start = math->iopack->get_comm();
  INIT_TIMER;
  START_TIMER;
  math->adder(dim, x, y, z, bw_x);
  STOP_TIMER("Adder");
  uint64_t comm_end = math->iopack->get_comm();
  double comm = (comm_end - comm_start) / (1.0 * (1ULL << 20)); // In MB
  cout << "Total Communication for Adder: " << comm << " mb" << endl;

  if (party == ALICE) {
    iopack->io->send_data(x, dim * sizeof(uint64_t));
    iopack->io->send_data(y, dim * sizeof(uint64_t));
    iopack->io->send_data(z, dim * sizeof(uint64_t));
  } else { // party == BOB
    uint64_t *x0 = new uint64_t[dim];
    uint64_t *y0 = new uint64_t[dim];
    uint64_t *z0 = new uint64_t[dim];
    iopack->io->recv_data(x0, dim * sizeof(uint64_t));
    iopack->io->recv_data(y0, dim * sizeof(uint64_t));
    iopack->io->recv_data(z0, dim * sizeof(uint64_t));

    for (int i = 0; i < dim; i++) {
      uint64_t X = (x[i] + x0[i]) & mask_x;
      uint64_t Y = (y[i] + y0[i]) & mask_x;
      int64_t D = signed_val(X - Y, bw_x);
      int64_t Z = signed_val(z[i] + z0[i], bw_x);
      int64_t expectedZ = -D;
      if (D < 0)
        expectedZ = D;
      assert(Z == expectedZ);
    }
    cout << "Adder Tests Passed" << endl;

    delete[] x0;
    delete[] y0;
    delete[] z0;
  }
  delete[] z;
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
  uint64_t *y = new uint64_t[dim];

  prg.random_data(x, dim * sizeof(uint64_t));
  prg.random_data(y, dim * sizeof(uint64_t));

  for (int i = 0; i < dim; i++) {
    x[i] &= mask_x;
    y[i] &= mask_x;
  }

  test_adder(x, y);

  delete[] x;
  delete[] y;
  delete math;
}
