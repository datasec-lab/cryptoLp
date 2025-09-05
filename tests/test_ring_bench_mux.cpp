#include "BuildingBlocks/aux-protocols.h"
#include "utils/emp-tool.h"
#include <iostream>

using namespace sci;
using namespace std;

int party, port = 32000;
string address = "127.0.0.1";
IOPack *iopack;
OTPack *otpack;
AuxProtocols *aux;

int dim = 1ULL << 16;
int bw_x = 32;

uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

void test_mux(uint8_t *sel, uint64_t *x, uint64_t *y) {
  uint64_t *z = new uint64_t[dim];

  uint64_t comm_start = aux->iopack->get_comm();
  INIT_TIMER;
  START_TIMER;
  aux->MUX(sel, x, y, z, dim, bw_x, bw_x);
  STOP_TIMER("MUX:");
  uint64_t comm_end = aux->iopack->get_comm();
  double comm = (comm_end - comm_start) / (1.0 * (1ULL << 20));
  cout << "Total Communication for MUX: " << comm << " mb" << endl;

  if (party == ALICE) {
    iopack->io->send_data(sel, dim * sizeof(uint8_t));
    iopack->io->send_data(x, dim * sizeof(uint64_t));
    iopack->io->send_data(y, dim * sizeof(uint64_t));
    iopack->io->send_data(z, dim * sizeof(uint64_t));
  } else { // party == BOB
    uint8_t *sel0 = new uint8_t[dim];
    uint64_t *x0 = new uint64_t[dim];
    uint64_t *y0 = new uint64_t[dim];
    uint64_t *z0 = new uint64_t[dim];
    iopack->io->recv_data(sel0, dim * sizeof(uint8_t));
    iopack->io->recv_data(x0, dim * sizeof(uint64_t));
    iopack->io->recv_data(y0, dim * sizeof(uint64_t));
    iopack->io->recv_data(z0, dim * sizeof(uint64_t));

    for (int i = 0; i < dim; i++) {
      uint8_t SEL = sel[i] ^ sel0[i];
      uint64_t X = (x[i] + x0[i]) & mask_x;
      uint64_t Y = (y[i] + y0[i]) & mask_x;
      uint64_t Z = (z[i] + z0[i]) & mask_x;
      uint64_t EXPECTED_Z = SEL ? X : Y;
      assert(Z == EXPECTED_Z);
    }

    cout << "MUX Tests Passed" << endl;

    delete[] sel0;
    delete[] x0;
    delete[] y0;
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

  aux = new AuxProtocols(party, iopack, otpack);

  PRG128 prg;

  uint8_t *sel = new uint8_t[dim]; // selector
  uint64_t *x = new uint64_t[dim];
  uint64_t *y = new uint64_t[dim];

  prg.random_data(sel, dim * sizeof(uint8_t));
  prg.random_data(x, dim * sizeof(uint64_t));
  prg.random_data(y, dim * sizeof(uint64_t));

  for (int i = 0; i < dim; i++) {
    sel[i] = sel[i] & 1;
    x[i] &= mask_x;
    y[i] &= mask_x;
  }

  test_mux(sel, x, y);

  delete[] sel;
  delete[] x;
  delete[] y;
  delete aux;
}
