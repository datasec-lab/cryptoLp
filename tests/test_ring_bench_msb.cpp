#include "BuildingBlocks/aux-protocols.h"
#include "utils/emp-tool.h"
#include <iostream>
using namespace sci;
using namespace std;

int party, port = 8000, dim = 1ULL << 16;
string address = "127.0.0.1";
IOPack *iopack;
OTPack *otpack;
AuxProtocols *aux;

int bw_x = 32;
uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

void test_msb(uint64_t *x) {
  int32_t shift = bw_x - 1;

  uint64_t shift_mask = (shift == 64 ? -1 : ((1ULL << shift) - 1));

  uint8_t *y = new uint8_t[dim];

  uint64_t comm_start = iopack->get_comm();
  INIT_TIMER;
  START_TIMER;
  aux->MSB(x, y, dim, bw_x);
  STOP_TIMER("MSB");
  double total_comm = (iopack->get_comm() - comm_start) / (1.0 * (1ULL << 20)); // In MB
  cout << "Total Communication for MSB: " << total_comm << " mb" << endl;
  
  if (party == ALICE) {
    iopack->io->send_data(x, dim * sizeof(uint64_t));
    iopack->io->send_data(y, dim * sizeof(uint8_t)); 
  } else {
    uint64_t *x0 = new uint64_t[dim];
    uint8_t *y0 = new uint8_t[dim];
    iopack->io->recv_data(x0, dim * sizeof(uint64_t));
    iopack->io->recv_data(y0, dim * sizeof(uint8_t));

    for (int i = 0; i < dim; i++) {
      int64_t X = signed_val(x[i] + x0[i], bw_x);
      uint8_t msb = (X >> shift) & 1;
      uint8_t Y = int(y[i] ^ y0[i]);
      assert(Y == msb);
    }
    cout << "MSB Tests passed" << endl;

    delete[] x0;
    delete[] y0;
  }
  delete[] y;
}

int main(int argc, char **argv) {
  ArgMapping amap;
  amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
  amap.arg("p", port, "Port Number");
  amap.arg("d", dim, "Size of vector");
  amap.arg("ip", address, "IP Address of server (ALICE)");
  amap.parse(argc, argv);

  iopack = new IOPack(party, port, "127.0.0.1");
  otpack = new OTPack(iopack, party);

  aux = new AuxProtocols(party, iopack, otpack);

  PRG128 prg;
  uint64_t *x = new uint64_t[dim];

  prg.random_data(x, dim * sizeof(uint64_t));
  for (int i = 0; i < dim; i++) {
    x[i] = x[i] & mask_x;
  }

  test_msb(x);

  delete[] x;
  return 0;
}
