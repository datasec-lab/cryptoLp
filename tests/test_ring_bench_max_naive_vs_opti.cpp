#include "Math/math-functions.h"
#include "utils/emp-tool.h"
#include <iostream>
#include <fstream>

using namespace sci;
using namespace std;

int party, port = 32000;
string address = "127.0.0.1";
IOPack *iopack;
OTPack *otpack;
MathFunctions *math;

#define NAIVE 0
#define OPTI 1

int dim;
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

void write_to_csv(double comm, double time, int mode) {
    std::ofstream file;
    if (mode == NAIVE) {
        if (party == ALICE) {
            file.open("max_naive_alice.csv", std::ios_base::app);
        } else {
            file.open("max_naive_bob.csv", std::ios_base::app);
        }
    } else {
        if (party == ALICE) {
            file.open("max_opti_alice.csv", std::ios_base::app);
        } else {
            file.open("max_opti_bob.csv", std::ios_base::app);
        }
    }
    if (file.is_open()) {
        file << comm << "," << time << std::endl;
        file.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
}

void test_max(uint64_t *x, int mode) {
  uint64_t *y = new uint64_t[1];

  uint64_t comm_start = math->iopack->get_comm();
  INIT_TIMER;
  START_TIMER;
  if (mode == NAIVE)
    math->max_naive(dim, x, y, bw_x);
  else
    math->max(dim, x, y, bw_x);
  auto end_timer = std::chrono::high_resolution_clock::now();
  uint64_t elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_timer - start_timer).count() + pause_timer;
  uint64_t comm_end = math->iopack->get_comm();
  double comm = (comm_end - comm_start) / (1.0 * (1ULL << 20));

  // Write to CSV
  write_to_csv(comm, elapsed_time, mode);

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

  cout << "Testing Max Naive" << endl;
  for (int i = 1; i <= 16; i++) {
    dim = 1 << i;

    uint64_t *x = new uint64_t[dim];

    prg.random_data(x, dim * sizeof(uint64_t));

    for (int i = 0; i < dim; i++) {
      x[i] &= mask_x;
    }
    math->ABS(dim, x, x, bw_x);

    test_max(x, NAIVE);

    delete[] x;
  }

  cout << "<><><><><><><><><><><><><><><><><><><><><><><><>" << endl;

  cout << "Testing Max Optimized" << endl;
  for (int i = 1; i <= 16; i++) {
    dim = 1 << i;

    uint64_t *x = new uint64_t[dim];

    prg.random_data(x, dim * sizeof(uint64_t));

    for (int i = 0; i < dim; i++) {
      x[i] &= mask_x;
    }
    math->ABS(dim, x, x, bw_x);

    test_max(x, OPTI);
  }

  delete math;
}

