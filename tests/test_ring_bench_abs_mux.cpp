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

void test_abs_mux() {
  int bw_x = 4, bw_y = 4;
  PRG128 prg;
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));

  uint8_t *sel = new uint8_t[dim]; // selector
  uint64_t *x = new uint64_t[dim];
  uint64_t *y = new uint64_t[dim];

  prg.random_data(sel, dim * sizeof(uint8_t));
  prg.random_data(x, dim * sizeof(uint64_t));
  for (int i = 0; i < dim; i++) {
    sel[i] = sel[i] & 1;
    x[i] = x[i] & mask_x;
  }

  uint64_t comm_start = iopack->get_comm();
  INIT_TIMER;
  START_TIMER;
  aux->ABS_multiplexer(sel, x, y, dim, bw_x, bw_y);
  STOP_TIMER("ABS Multiplexer");
  double total_comm = (iopack->get_comm() - comm_start) / (1.0 * (1ULL << 20)); // In MB
  cout << "Total Communication for ABS Multiplexer: " << total_comm << " mb" << endl;
  
  if (party == ALICE) {
    iopack->io->send_data(sel, dim * sizeof(uint8_t));
    iopack->io->send_data(x, dim * sizeof(uint64_t));
    iopack->io->send_data(y, dim * sizeof(uint64_t)); 
  } else {
    uint8_t *sel0 = new uint8_t[dim];
    uint64_t *x0 = new uint64_t[dim];
    uint64_t *y0 = new uint64_t[dim];
    iopack->io->recv_data(sel0, dim * sizeof(uint8_t));
    iopack->io->recv_data(x0, dim * sizeof(uint64_t));
    iopack->io->recv_data(y0, dim * sizeof(uint64_t));

    for (int i = 0; i < dim; i++) {
      assert((((1 - 2*uint64_t(sel0[i] ^ sel[i])) * (x0[i] + x[i])) & mask_y) ==
             ((y0[i] + y[i]) & mask_y));
    }
    cout << "ABS MUX Tests passed" << endl;

    delete[] sel0;
    delete[] x0;
    delete[] y0;
  }
  delete[] sel;
  delete[] x;
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

  test_abs_mux();

  return 0;
}
