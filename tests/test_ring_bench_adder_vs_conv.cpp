#include "Math/math-functions.h"
#include "utils/emp-tool.h"
#include <iostream>
#include <fstream>
#include <chrono>

using namespace sci;
using namespace std;

int party, port = 32000;
string address = "127.0.0.1";
IOPack *iopack;
OTPack *otpack;
MathFunctions *math;

#define ADDER 1
#define CONV 0

int dim;
int bw_x = 32;

uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

void write_to_csv(double comm, double time, int mode) {
    std::ofstream file;
    if (mode == ADDER) {
        if (party == ALICE) {
            file.open("adder_alice.csv", std::ios_base::app);
        } else {
            file.open("adder_bob.csv", std::ios_base::app);
        }
    } else {
        if (party == ALICE) {
            file.open("conv_alice.csv", std::ios_base::app);
        } else {
            file.open("conv_bob.csv", std::ios_base::app);
        }
    }
    if (file.is_open()) {
        file << comm << "," << time << std::endl;
        file.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
}

std::pair<double, double> evaluate_performance(uint64_t *x, uint64_t *y, int mode) {
  uint64_t *z = new uint64_t[dim];

  uint64_t comm_start = math->iopack->get_comm();
  INIT_TIMER;
  START_TIMER;
  if (mode == ADDER)
    math->adder(dim, x, y, z, bw_x);
  else
    math->conv(dim, x, y, z, bw_x, bw_x, bw_x);
  auto end_timer = std::chrono::high_resolution_clock::now();
  uint64_t elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_timer - start_timer).count() + pause_timer;
  uint64_t comm_end = math->iopack->get_comm();
  double comm = (comm_end - comm_start) / (1.0 * (1ULL << 20));

  delete[] z;

  return std::make_pair(comm, elapsed_time);
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

  const int num_runs = 10;
  PRG128 prg;

  cout << "Testing Adder" << endl;
  for (int i = 0; i <= 16; i++) {
    dim = 1 << i;

    vector<double> comms, times;

    for (int run = 0; run < num_runs; run++) {
      uint64_t *x = new uint64_t[dim];
      uint64_t *y = new uint64_t[dim];

      prg.random_data(x, dim * sizeof(uint64_t));

      for (int i = 0; i < dim; i++) {
        x[i] &= mask_x;
        y[i] &= mask_x;
      }

      auto result = evaluate_performance(x, y, ADDER);

      comms.push_back(result.first);
      times.push_back(result.second);

      delete[] x;
      delete[] y;
    }

    double avg_comm = accumulate(comms.begin(), comms.end(), 0.0) / num_runs;
    double avg_time = accumulate(times.begin(), times.end(), 0.0) / num_runs;

    write_to_csv(avg_comm, avg_time, ADDER);
  }

  cout << "<><><><><><><><><><><><><><><><><><><><><><><><>" << endl;

  cout << "Testing Conv" << endl;
  for (int i = 0; i <= 16; i++) {
    dim = 1 << i;

    vector<double> comms, times;

    for (int run = 0; run < num_runs; run++) {
      uint64_t *x = new uint64_t[dim];
      uint64_t *y = new uint64_t[dim];

      prg.random_data(x, dim * sizeof(uint64_t));

      for (int i = 0; i < dim; i++) {
        x[i] &= mask_x;
        y[i] &= mask_x;
      }

      auto result = evaluate_performance(x, y, CONV);

      comms.push_back(result.first);
      times.push_back(result.second);

      delete[] x;
      delete[] y;
    }

    double avg_comm = accumulate(comms.begin(), comms.end(), 0.0) / num_runs;
    double avg_time = accumulate(times.begin(), times.end(), 0.0) / num_runs;

    write_to_csv(avg_comm, avg_time, CONV);
  }

  delete math;
}

