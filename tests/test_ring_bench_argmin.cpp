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

int dim = 1ULL << 16;
int bw_x = 32;

uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

std::pair<uint64_t, uint64_t> find_min(uint64_t *x, int dim) {
  if (dim == 0) {
    throw std::invalid_argument("find_min: dim must be greater than 0");
  }

  int n = dim;
  uint64_t m;
  uint64_t m_index;
  while (n > 1) {
    int new_dim = (n + 1) / 2;
    uint64_t *min = new uint64_t[new_dim];
    uint64_t *min_index = new uint64_t[new_dim];
    for (int i = 0; i < new_dim; i++) {
      if (!(n % 2)) {
        if (x[2 * i] < x[2 * i + 1]) {
          min[i] = x[2 * i], min_index[i] = 2 * i;
        } else {
          min[i] = x[2 * i + 1], min_index[i] = 2 * i + 1;
        }
      } else {
        if (i != new_dim - 1) {
          if (x[2 * i] < x[2 * i + 1]) {
            min[i] = x[2 * i], min_index[i] = 2 * i;
          } else {
            min[i] = x[2 * i + 1], min_index[i] = 2 * i + 1;
          }
        } else {
          min[i] = x[2 * i], min_index[i] = 2 * i;
        }
      }
    }
    n = new_dim;
    if (n == 1) {
      m = min[0];
      m_index = min_index[0];
    }
  }
  return std::make_pair(m, m_index);
}

void test_argmin(uint64_t *x) {
  uint64_t *y = new uint64_t[1];
  uint64_t *yi = new uint64_t[1];

  uint64_t comm_start = math->iopack->get_comm();
  INIT_TIMER;
  START_TIMER;
  math->argmin(dim, x, y, yi, bw_x);

  STOP_TIMER("ArgMin:");
  uint64_t comm_end = math->iopack->get_comm();
  double comm = (comm_end - comm_start) / (1.0 * (1ULL << 20));
  cout << "Total Communication for ArgMin: " << comm << " mb" << endl;

  if (party == ALICE) {
    iopack->io->send_data(x, dim * sizeof(uint64_t));
    iopack->io->send_data(y, sizeof(uint64_t));
    iopack->io->send_data(yi, sizeof(uint64_t));
  } else { // party == BOB
    uint64_t *x0 = new uint64_t[dim];
    uint64_t *y0 = new uint64_t[1];
    uint64_t *yi0 = new uint64_t[1];
    iopack->io->recv_data(x0, dim * sizeof(uint64_t));
    iopack->io->recv_data(y0, sizeof(uint64_t));
    iopack->io->recv_data(yi0, sizeof(uint64_t));

    uint64_t *temp_x = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
      temp_x[i] = (x[i] + x0[i]) & mask_x;
      cout << "temp_x[" << i << "]: " << temp_x[i] << endl;
    }

    std::pair<uint64_t, uint64_t> result = find_min(temp_x, dim);
    uint64_t min = result.first;
    uint64_t min_index = result.second;
    
    uint64_t Y = (y[0] + y0[0]) & mask_x;
    uint64_t YI = (yi[0] + yi0[0]) & mask_x;

    // cout << "Y: " << Y << endl;
    cout << "YI: " << YI << endl;
    // cout << "min: " << min << endl;
    // cout << "min_index: " << min_index << endl;
    
    assert(Y == min);
    assert(YI == min_index);
    

    cout << "ArgMin Tests Passed" << endl;

    delete[] x0;
    delete[] y0;
    delete[] yi0;
    delete[] temp_x;
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

  test_argmin(x);

  delete[] x;
  delete math;
}
