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

int m = 2; // m
uint8_t mask_m = (1 << m) - 1;
int num_digits = ceil((double)bw_x / m); // q

void test_mill(uint64_t *x) {
  uint8_t *y = new uint8_t[dim];
  uint8_t *z = new uint8_t[(num_digits - 1) * dim];

  uint64_t comm_start = iopack->get_comm();
  INIT_TIMER;
  START_TIMER;
  // math->millionaire(x, y, dim, bw_x);
  // math->millionaire_ext_naive(x, y, dim, bw_x);
  // math->millionaire_ext(x, y, dim, bw_x);
  math->millionaire_ext_plus(x, y, dim, bw_x);
  STOP_TIMER("Mill");
  double total_comm = (iopack->get_comm() - comm_start) / (1.0 * (1ULL << 20)); // In MB
  cout << "Total Communication for Mill: " << total_comm << " mb" << endl;

  if (party == ALICE) {
    iopack->io->send_data(x, dim * sizeof(uint64_t));
    iopack->io->send_data(y, dim * sizeof(uint8_t));
    iopack->io->send_data(z, (num_digits - 1) * dim * sizeof(uint8_t));
  } else { // party == BOB
    uint64_t *x0 = new uint64_t[dim];
    uint8_t *y0 = new uint8_t[dim];
    uint8_t *z0 = new uint8_t[(num_digits - 1) * dim];
    iopack->io->recv_data(x0, dim * sizeof(uint64_t));
    iopack->io->recv_data(y0, dim * sizeof(uint8_t));
    iopack->io->recv_data(z0, (num_digits - 1) * dim * sizeof(uint8_t));

    uint8_t *res = new uint8_t[dim];
    uint8_t *lt = new uint8_t[num_digits * dim];
    uint8_t *eq = new uint8_t[num_digits * dim];
    uint8_t *digits_A = new uint8_t[num_digits * dim];
    uint8_t *digits_B = new uint8_t[num_digits * dim];

    for (int i = 0; i < dim; i++) {
      res[i] = (x0[i] < x[i]);

      for (int j = 0; j < num_digits; j++) {
        digits_A[i * num_digits + j] = (uint8_t)(x0[i] >> j * m) & mask_m;
        digits_B[i * num_digits + j] = (uint8_t)(x[i] >> j * m) & mask_m;
      }

      for (int j = 0; j < num_digits; j++) {
        lt[i * num_digits + j] = (digits_A[i * num_digits + j] < digits_B[i * num_digits + j]);
        eq[i * num_digits + j] = (digits_A[i * num_digits + j] == digits_B[i * num_digits + j]);
      }
      
      for (int j = 0; j < num_digits; j++) {
        if (j != num_digits - 1) {
          uint8_t Z = z0[i * (num_digits - 1) + j] ^ z[i * (num_digits - 1) + j];
          uint8_t lt_eq = lt[i * num_digits + j] & eq[i * num_digits + (j + 1)];
        }
      }

      // cout << "x0: " << x0[i] << endl;
      // cout << "x: " << x[i] << endl;

      // if (res[i] != y0[i] ^ y[i]) {
      //   cout << "x0: " << x0[i] << endl;
      //   cout << "x: " << x[i] << endl;
      // }
      assert(res[i] == y0[i] ^ y[i]);
    }
    cout << "Mill Tests Passed" << endl;

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

  test_mill(x);

  delete[] x;
  delete math;
}
