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
int bw_z = 32;
// int s_x = 12;
// int s_y = 12;
// int s_z = 24;

uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));
uint64_t mask_z = (bw_z == 64 ? -1 : ((1ULL << bw_z) - 1));

uint64_t computeULPErr(double calc, double actual, int SCALE) {
  int64_t calc_fixed = (double(calc) * (1ULL << SCALE));
  int64_t actual_fixed = (double(actual) * (1ULL << SCALE));
  uint64_t ulp_err = (calc_fixed - actual_fixed) > 0
                         ? (calc_fixed - actual_fixed)
                         : (actual_fixed - calc_fixed);
  return ulp_err;
}

void test_conv(uint64_t *x, uint64_t *y) {
  uint64_t *z = new uint64_t[dim];

  /* Test case
  * x -> 8.25
  * y -> 3.5
  * xy -> 28.875
  */
  // if (party == ALICE) {
  //   x[0] = 22;
  //   y[0] = 14;
  // } else { // party == BOB
  //   x[0] = 44;
  //   y[0] = 14;
  // }

  uint64_t comm_start = math->iopack->get_comm();
  INIT_TIMER;
  START_TIMER;
  math->conv(dim, x, y, z, bw_x, bw_y, bw_z);
  STOP_TIMER("Conv:");
  uint64_t comm_end = math->iopack->get_comm();
  double comm = (comm_end - comm_start) / (1.0 * (1ULL << 20));
  cout << "Total Communication for Conv: " << comm << " mb" << endl;

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

    uint64_t total_err = 0;
    uint64_t max_ULP_err = 0;
    for (int i = 0; i < dim; i++) {
      uint64_t X = (x0[i] + x[i]) & mask_x;
      uint64_t Y = (y0[i] + y[i]) & mask_y;
      uint64_t P = (X * Y) & mask_z;
      uint64_t Z = (z0[i] + z[i]) & mask_z;

      assert(Z == P);


      // // ############## Compute ULP Error ##############
      // double dbl_x = X / double(1ULL << s_x);
      // double dbl_y = Y / double(1ULL << s_y);
      // uint64_t prod = uint64_t((dbl_x * dbl_y) * (1ULL << s_z)) & mask_z;
      // double dbl_prod = prod / double(1ULL << s_z);
      // double dbl_z = Z / double(1ULL << s_z);
      // cout << "dbl_prod: " << dbl_prod << endl;
      // cout << "dbl_z: " << Z / double(1ULL << s_z) << endl;
      // uint64_t err = computeULPErr(dbl_z, dbl_prod, s_z);
      // total_err += err;
      // max_ULP_err = std::max(max_ULP_err, err);
    }

    cout << "Conv Tests Passed" << endl;
    // cerr << "Average ULP error: " << total_err / dim << endl;
    // cerr << "Max ULP error: " << max_ULP_err << endl;

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
    y[i] &= mask_y;
  }

  test_conv(x, y);

  delete[] x;
  delete[] y;
  delete math;
}
