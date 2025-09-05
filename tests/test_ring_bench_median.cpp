#include "Math/math-functions.h"
#include "utils/emp-tool.h"
#include <iostream>
#include <utility>

using namespace sci;
using namespace std;

int party, port = 32000;
string address = "127.0.0.1";
IOPack *iopack;
OTPack *otpack;
MathFunctions *math;

int dim = 1ULL << 4;
int bw_x = 32;

uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

// sort in ascending order
void sort(uint64_t *x, int dim) {
  int64_t *x_signed = new int64_t[dim];
  for (int i = 0; i < dim; i++) {
    x_signed[i] = signed_val(x[i], bw_x);
  }
  for (int i = 0; i < dim - 1; i++) {
    for (int j = i + 1; j < dim; j++) {
      if (x_signed[i] > x_signed[j]) {
        uint64_t temp = x[i];
        int64_t temp_signed = x_signed[i];
        x[i] = x[j];
        x_signed[i] = x_signed[j];
        x[j] = temp;
        x_signed[j] = temp_signed;
      }
    }
  }
  delete[] x_signed;
}

uint64_t median(uint64_t *x) {
  sort(x, dim);
  if (dim % 2 == 0) {
    int64_t tmp0 = signed_val(x[dim / 2 - 1], bw_x) / 2;
    int64_t tmp1 = signed_val(x[dim / 2], bw_x) / 2;
    return unsigned_val(tmp0 + tmp1, bw_x);
  } else {
    return x[dim / 2];
  }
}

void test_median(uint64_t *x) {
  uint64_t comm_start = math->iopack->get_comm();
  uint64_t *y = new uint64_t[1];
  
  INIT_TIMER;
  START_TIMER;
  math->median(dim, x, y, bw_x);
  STOP_TIMER("Median:");
  uint64_t comm_end = math->iopack->get_comm();
  double comm = (comm_end - comm_start) / (1.0 * (1ULL << 20));
  cout << "Total Communication for median: " << comm << " mb" << endl;

  if (party == ALICE) {
    iopack->io->send_data(x, dim * sizeof(uint64_t));
    iopack->io->send_data(y, 1 * sizeof(uint64_t));
  } else { // party == BOB
    uint64_t *x0 = new uint64_t[dim];
    uint64_t *y0 = new uint64_t[1];
    iopack->io->recv_data(x0, dim * sizeof(uint64_t));
    iopack->io->recv_data(y0, 1 * sizeof(uint64_t));

    for (int i = 0; i < dim; i++) {
      x[i] = (x[i] + x0[i]) & mask_x;
    }

    uint64_t median_expected = median(x);

    uint64_t median;
    median = (y[0] + y0[0]) & mask_x;
    
    assert(median == median_expected);

    cout << "Median Tests Passed" << endl;

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
  uint64_t* *temp;
  temp = new uint64_t*[dim];

  cout << "Preparing test data..." << endl;
  for (int i = 0; i < dim; i++) {
    if (party == ALICE) {
      temp[i] = new uint64_t[1];
      prg.random_data(temp[i], sizeof(uint64_t));
      temp[i][0] &= mask_x;
      x[i] = temp[i][0];
      iopack->io->send_data(temp[i], sizeof(uint64_t));
    } else {
      uint64_t *temp0 = new uint64_t[1];
      iopack->io->recv_data(temp0, sizeof(uint64_t));
      bool flag = true;
      while (flag) {
        temp[i] = new uint64_t[1];
        prg.random_data(temp[i], sizeof(uint64_t));
        temp[i][0] &= mask_x;
        int64_t temp_signed = signed_val(temp[i][0] + temp0[0], bw_x);
        if (temp_signed >= -(1 << (bw_x -2)) && temp_signed <= (1 << (bw_x - 2))) {
          flag = false;
          x[i] = temp[i][0];
        }
      }
    }
  }
  cout << "Done" << endl;

  for (int i = 0; i < dim; i++) {
    delete[] temp[i];
  }

  test_median(x);

  delete[] x;
  delete[] temp;
  delete math;
}
