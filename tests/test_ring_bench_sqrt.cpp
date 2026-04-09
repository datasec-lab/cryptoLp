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
int bw_y = 32;
int s_x = 12;
int s_y = 12;

uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << (bw_x - 1)) - 1));

void test_sqrt(uint64_t *x) {
  uint64_t *y = new uint64_t[dim];

  uint64_t comm_start = math->iopack->get_comm();
  INIT_TIMER;
  START_TIMER;
  math->sqrt(dim, x, y, bw_x, bw_y, s_x, s_y, false);
  STOP_TIMER("Sqrt:");
  uint64_t comm_end = math->iopack->get_comm();
  double comm = (comm_end - comm_start) / (1.0 * (1ULL << 20));
  cout << "Total Communication for Sqrt: " << comm << " MB" << endl;

  if (party == ALICE) {
    iopack->io->send_data(y, dim * sizeof(uint64_t));
  } else {
    uint64_t *y0 = new uint64_t[dim];
    iopack->io->recv_data(y0, dim * sizeof(uint64_t));
    cout << "Sqrt Tests Passed" << endl;
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
  uint64_t comm_before = iopack->get_comm();
  INIT_TIMER;
  START_TIMER;
  otpack = new OTPack(iopack, party);

  STOP_TIMER("OTPack setup");
  uint64_t comm_after = iopack->get_comm();
  cout << "OTPack setup communication: "
       << (comm_after - comm_before) / (1024.0 * 1024.0) << " MB" << endl;
  math = new MathFunctions(party, iopack, otpack);

  PRG128 prg;
  uint64_t *x = new uint64_t[dim];

  if (party == ALICE) {
    prg.random_data(x, dim * sizeof(uint64_t));
    for (int i = 0; i < dim; i++) {
      x[i] &= mask_x;
    }
  } else {
    for (int i = 0; i < dim; i++) {
      x[i] = 0;
    }
  }

  test_sqrt(x);

  delete[] x;
  delete math;
}