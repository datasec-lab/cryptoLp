/*
Authors: Deevashwer Rathee
Copyright:
Copyright (c) 2021 Microsoft Research
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "NonLinear/relu-ring.h"
#include "Math/math-functions.h"
#include "Millionaire/bit-triple-generator.h"

using namespace std;
using namespace sci;

#define KKOT_LIMIT 8
#define SQRT_LOOKUP_SCALE 2
#define RING 0

MathFunctions::MathFunctions(int party, IOPack *iopack, OTPack *otpack) {
  this->party = party;
  this->iopack = iopack;
  this->otpack = otpack;
  this->aux = new AuxProtocols(party, iopack, otpack);
  this->xt = new XTProtocol(party, iopack, otpack);
  this->trunc = new Truncation(party, iopack, otpack);
  this->mult = new LinearOT(party, iopack, otpack);
}

MathFunctions::~MathFunctions() {
  delete aux;
  delete xt;
  delete trunc;
  delete mult;
}

// A0 \in (1/4, 1)
uint64_t lookup_A0(uint64_t index, int m) {
  uint64_t k = 1ULL << m;
  double p = 1 + (double(index) / double(k));
  double A1 = 1.0 / (p * (p + 1.0 / double(k)));
  int32_t scale = m + 3;
  uint64_t mask = (1ULL << scale) - 1;
  uint64_t val = uint64_t(A1 * (1ULL << scale)) & mask;
  return val;
}

// A1 \in (1/2, 1)
uint64_t lookup_A1(uint64_t index, int m) {
  uint64_t k = 1ULL << m;
  double p = 1 + (double(index) / double(k));
  double z = (p * (p + (1.0 / double(k))));
  double A1 = ((1.0 / double(k * 2)) + sqrt(z)) / z;
  int32_t scale = 2 * m + 2;
  uint64_t mask = (1ULL << scale) - 1;
  uint64_t val = uint64_t(A1 * (1ULL << scale)) & mask;
  return val;
}

void MathFunctions::reciprocal_approximation(int32_t dim, int32_t m,
                                             uint64_t *dn, uint64_t *out,
                                             int32_t bw_dn, int32_t bw_out,
                                             int32_t s_dn, int32_t s_out) {
  assert(bw_out == m + s_dn + 4);
  assert(s_out == m + s_dn + 4);

  uint64_t s_dn_mask = (1ULL << s_dn) - 1;
  uint64_t m_mask = (1ULL << m) - 1;
  uint64_t s_min_m_mask = (1ULL << (s_dn - m)) - 1;

  uint64_t *tmp_1 = new uint64_t[dim];
  uint64_t *tmp_2 = new uint64_t[dim];

  for (int i = 0; i < dim; i++) {
    tmp_1[i] = dn[i] & s_dn_mask;
  }
  trunc->truncate_and_reduce(dim, tmp_1, tmp_2, s_dn - m, s_dn);

  int M = (1ULL << m);
  uint64_t c0_mask = (1ULL << (m + 4)) - 1;
  uint64_t c1_mask = (1ULL << (2 * m + 3)) - 1;
  uint64_t *c0 = new uint64_t[dim];
  uint64_t *c1 = new uint64_t[dim];
  if (party == ALICE) {
    uint64_t **spec;
    spec = new uint64_t *[dim];
    PRG128 prg;
    prg.random_data(c0, dim * sizeof(uint64_t));
    prg.random_data(c1, dim * sizeof(uint64_t));
    for (int i = 0; i < dim; i++) {
      spec[i] = new uint64_t[M];
      c0[i] &= c0_mask;
      c1[i] &= c1_mask;
      for (int j = 0; j < M; j++) {
        int idx = (tmp_2[i] + j) & m_mask;
        spec[i][j] = (lookup_A0(idx, m) - c0[i]) & c0_mask;
        spec[i][j] <<= (2 * m + 3);
        spec[i][j] |= (lookup_A1(idx, m) - c1[i]) & c1_mask;
      }
    }
    aux->lookup_table<uint64_t>(spec, nullptr, nullptr, dim, m, 3 * m + 7);

    for (int i = 0; i < dim; i++)
      delete[] spec[i];
    delete[] spec;
  } else {
    aux->lookup_table<uint64_t>(nullptr, tmp_2, c1, dim, m, 3 * m + 7);

    for (int i = 0; i < dim; i++) {
      c0[i] = (c1[i] >> (2 * m + 3)) & c0_mask;
      c1[i] = c1[i] & c1_mask;
    }
  }
  for (int i = 0; i < dim; i++) {
    tmp_1[i] = dn[i] & s_min_m_mask;
  }
  uint8_t *zero_shares = new uint8_t[dim];
  for (int i = 0; i < dim; i++) {
    zero_shares[i] = 0;
  }

  // Unsigned mult
  mult->hadamard_product(dim, c0, tmp_1, tmp_2, m + 4, s_dn - m, s_dn + 4,
                         false, false, MultMode::None, zero_shares, nullptr);

  xt->z_extend(dim, tmp_2, tmp_1, s_dn + 4, s_dn + m + 4, zero_shares);

  uint64_t out_mask = (1ULL << (s_dn + m + 4)) - 1;
  uint64_t scale_up = (1ULL << (s_dn - m + 1));
  for (int i = 0; i < dim; i++) {
    out[i] = ((c1[i] * scale_up) - tmp_1[i]) & out_mask;
  }

  delete[] tmp_1;
  delete[] tmp_2;
  delete[] c0;
  delete[] c1;
  delete[] zero_shares;
}

void MathFunctions::div(int32_t dim, uint64_t *nm, uint64_t *dn, uint64_t *out,
                        int32_t bw_nm, int32_t bw_dn, int32_t bw_out,
                        int32_t s_nm, int32_t s_dn, int32_t s_out,
                        bool signed_nm, bool compute_msnzb) {
  assert(s_out <= s_dn);

  // out_precision = iters * (2*m + 2)
  int32_t m, iters;
  m = (s_out <= 18 ? ceil((s_out - 2) / 2.0)
                   : ceil((ceil(s_out / 2.0) - 2) / 2.0));
  iters = (s_out <= 18 ? 0 : 1);

  int32_t s_tmp_dn;
  int32_t bw_adjust;
  int32_t s_adjust;
  uint64_t *tmp_dn;
  uint64_t *adjust;
  if (compute_msnzb) {
    s_tmp_dn = bw_dn - 1;
    bw_adjust = bw_dn + 1;
    s_adjust = bw_dn - 1 - s_dn;
    uint64_t mask_adjust = (bw_adjust == 64 ? -1 : ((1ULL << bw_adjust) - 1));
    // MSB is always 0, thus, not including it
    uint8_t *msnzb_vector_bool = new uint8_t[dim * bw_dn];
    uint64_t *msnzb_vector = new uint64_t[dim * bw_dn];
    aux->msnzb_one_hot(dn, msnzb_vector_bool, bw_dn, dim);
    aux->B2A(msnzb_vector_bool, msnzb_vector, dim * bw_dn, bw_adjust);
    // adjust: bw = bw_dn, scale = bw_dn - 1 - s_dn
    adjust = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
      adjust[i] = 0;
      for (int j = 0; j < bw_dn; j++) {
        adjust[i] += (1ULL << (bw_dn - 1 - j)) * msnzb_vector[i * bw_dn + j];
      }
      adjust[i] &= mask_adjust;
    }
    // tmp_dn: bw = bw_dn, scale = bw_dn - 1
    tmp_dn = new uint64_t[dim];
    mult->hadamard_product(dim, dn, adjust, tmp_dn, bw_dn, bw_dn + 1, bw_dn + 1,
                           false, false, MultMode::None);

    delete[] msnzb_vector_bool;
    delete[] msnzb_vector;
  } else {
    // tmp_dn: bw = s_dn + 1, scale = s_dn
    s_tmp_dn = s_dn;
    tmp_dn = dn;
  }

  uint64_t *tmp_1 = new uint64_t[dim];
  uint64_t *tmp_2 = new uint64_t[dim];
  // tmp_1: bw = s_tmp_dn + m + 4, scale = s_tmp_dn + m + 3
  reciprocal_approximation(dim, m, tmp_dn, tmp_1, bw_dn, s_tmp_dn + m + 4,
                           s_tmp_dn, s_tmp_dn + m + 4);

  uint64_t *w0 = new uint64_t[dim];
  // w0: bw = s_out + 1, scale = s_out
  trunc->truncate_and_reduce(dim, tmp_1, w0, s_tmp_dn + m + 3 - s_out,
                             s_tmp_dn + m + 4);

  uint8_t *msb_nm = new uint8_t[dim];
  aux->MSB(nm, msb_nm, dim, bw_nm);
  uint8_t *zero_shares = new uint8_t[dim];
  for (int i = 0; i < dim; i++) {
    zero_shares[i] = 0;
  }

  // a0: bw = bw_out, scale = s_out
  uint64_t *a0 = new uint64_t[dim];
  // Mixed mult with w0 unsigned
  mult->hadamard_product(dim, nm, w0, tmp_1, bw_nm, s_out + 1, s_out + bw_nm,
                         signed_nm, false, MultMode::None, msb_nm, zero_shares);
  trunc->truncate_and_reduce(dim, tmp_1, tmp_2, s_nm, s_out + bw_nm);
  if ((bw_nm - s_nm) >= (bw_out - s_out)) {
    aux->reduce(dim, tmp_2, a0, bw_nm - s_nm + s_out, bw_out);
  } else {
    if (signed_nm) {
      xt->s_extend(dim, tmp_2, a0, s_out + bw_nm - s_nm, bw_out, msb_nm);
    } else {
      xt->z_extend(dim, tmp_2, a0, s_out + bw_nm - s_nm, bw_out, nullptr);
    }
  }

  if (compute_msnzb) {
    int32_t bw_tmp1 =
        (bw_out + s_adjust < bw_adjust ? bw_adjust : bw_out + s_adjust);
    // tmp_1: bw = bw_tmp1, scale = s_out + s_adjust
    mult->hadamard_product(dim, a0, adjust, tmp_1, bw_out, bw_adjust, bw_tmp1,
                           signed_nm, false, MultMode::None,
                           (signed_nm ? msb_nm : nullptr), zero_shares);
    // a0: bw = bw_out, scale = s_out
    trunc->truncate_and_reduce(dim, tmp_1, a0, s_adjust, bw_out + s_adjust);
  }

  // tmp_1: bw = s_tmp_dn + 2, scale = s_tmp_dn
  uint64_t s_plus_2_mask = (1ULL << (s_tmp_dn + 2)) - 1;
  for (int i = 0; i < dim; i++) {
    tmp_1[i] = tmp_dn[i] & s_plus_2_mask;
  }

  if (iters > 0) {
    // d0: bw = s_out + 2, scale = s_out
    uint64_t *d0 = new uint64_t[dim];
    mult->hadamard_product(dim, w0, tmp_1, tmp_2, s_out + 1, s_tmp_dn + 2,
                           s_out + s_tmp_dn + 2, false, false, MultMode::None,
                           zero_shares, zero_shares);
    trunc->truncate_and_reduce(dim, tmp_2, d0, s_tmp_dn, s_out + s_tmp_dn + 2);

    // e0: bw = s_out + 2, scale = s_out
    // e0 = 1 - d0
    uint64_t *e0 = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
      e0[i] = (party == ALICE ? (1ULL << (s_out)) : 0) - d0[i];
    }

    uint64_t e_mask = (1ULL << (s_out + 2)) - 1;
    uint64_t *a_curr = new uint64_t[dim];
    uint64_t *e_curr = new uint64_t[dim];
    uint64_t *a_prev = a0;
    uint64_t *e_prev = e0;
    for (int i = 0; i < iters - 1; i++) {
      // tmp_1: bw = 2*s_out+2, scale: 2*s_out
      mult->hadamard_product(dim, e_prev, e_prev, tmp_1, s_out + 2, s_out + 2,
                             2 * s_out + 2, true, true, MultMode::None,
                             zero_shares, zero_shares);
      // e_curr: bw = s_out + 2, scale: s_out
      trunc->truncate_and_reduce(dim, tmp_1, e_curr, s_out, 2 * s_out + 2);
      // e_prev = 1 + e_prev
      for (int j = 0; j < dim; j++) {
        e_prev[j] =
            ((party == ALICE ? (1ULL << (s_out)) : 0) + e_prev[j]) & e_mask;
      }
      // tmp_1: bw = bw_out + s_out, scale: 2*s_out
      mult->hadamard_product(dim, a_prev, e_prev, tmp_1, bw_out, s_out + 2,
                             bw_out + s_out, signed_nm, true, MultMode::None,
                             (signed_nm ? msb_nm : nullptr), zero_shares);
      // a_curr: bw = bw_out, scale: s_out
      trunc->truncate_and_reduce(dim, tmp_1, a_curr, s_out, bw_out + s_out);
      memcpy(a_prev, a_curr, dim * sizeof(uint64_t));
      memcpy(e_prev, e_curr, dim * sizeof(uint64_t));
    }
    // e_prev = 1 + e_prev
    for (int j = 0; j < dim; j++) {
      e_prev[j] =
          ((party == ALICE ? (1ULL << (s_out)) : 0) + e_prev[j]) & e_mask;
    }
    // out: bw = bw_out, scale: s_out
    // Mixed mult with e_prev unsigned
    mult->hadamard_product(dim, a_prev, e_prev, tmp_1, bw_out, s_out + 2,
                           bw_out + s_out, signed_nm, false, MultMode::None,
                           (signed_nm ? msb_nm : nullptr), zero_shares);
    trunc->truncate_and_reduce(dim, tmp_1, out, s_out, bw_out + s_out);

    delete[] d0;
    delete[] e0;
    delete[] a_curr;
    delete[] e_curr;
  } else {
    memcpy(out, a0, dim * sizeof(uint64_t));
  }

  delete[] tmp_1;
  delete[] tmp_2;
  delete[] w0;
  delete[] a0;
  delete[] msb_nm;
  delete[] zero_shares;
  if (compute_msnzb) {
    delete[] tmp_dn;
    delete[] adjust;
  }
}

uint64_t lookup_neg_exp(uint64_t val_in, int32_t s_in, int32_t s_out) {
  if (s_in < 0) {
    s_in *= -1;
    val_in *= (1ULL << (s_in));
    s_in = 0;
  }
  uint64_t res_val =
      exp(-1.0 * (val_in / double(1ULL << s_in))) * (1ULL << s_out);
  return res_val;
}

void MathFunctions::lookup_table_exp(int32_t dim, uint64_t *x, uint64_t *y,
                                     int32_t bw_x, int32_t bw_y, int32_t s_x,
                                     int32_t s_y) {
  assert(bw_y >= (s_y + 2));

  int LUT_size = KKOT_LIMIT;

  uint64_t bw_x_mask = (bw_x == 64 ? -1 : (1ULL << bw_x) - 1);
  uint64_t LUT_out_mask = ((s_y + 2) == 64 ? -1 : (1ULL << (s_y + 2)) - 1);

  uint64_t *tmp_1 = new uint64_t[dim];
  for (int i = 0; i < dim; i++) {
    tmp_1[i] = (-1 * x[i]) & bw_x_mask;
  }
  int digit_size = LUT_size;
  int num_digits = ceil(double(bw_x) / digit_size);
  int last_digit_size = bw_x - (num_digits - 1) * digit_size;
  uint64_t *x_digits = new uint64_t[num_digits * dim];

  aux->digit_decomposition_sci(dim, tmp_1, x_digits, bw_x, digit_size);

  uint64_t digit_mask = (digit_size == 64 ? -1 : (1ULL << digit_size) - 1);
  uint64_t last_digit_mask =
      (last_digit_size == 64 ? -1 : (1ULL << last_digit_size) - 1);
  int N = (1ULL << digit_size);
  int last_N = (1ULL << last_digit_size);
  int N_digits = (digit_size == last_digit_size ? num_digits : num_digits - 1);
  uint64_t *digits_exp = new uint64_t[num_digits * dim];
  if (party == ALICE) {
    uint64_t **spec;
    spec = new uint64_t *[num_digits * dim];
    PRG128 prg;
    prg.random_data(digits_exp, num_digits * dim * sizeof(uint64_t));
    for (int i = 0; i < N_digits * dim; i++) {
      int digit_idx = i / dim;
      spec[i] = new uint64_t[N];
      digits_exp[i] &= LUT_out_mask;
      for (int j = 0; j < N; j++) {
        int idx = (x_digits[i] + j) & digit_mask;
        spec[i][j] = (lookup_neg_exp(idx, s_x - digit_size * digit_idx, s_y) -
                      digits_exp[i]) &
                     LUT_out_mask;
      }
    }
    aux->lookup_table<uint64_t>(spec, nullptr, nullptr, N_digits * dim,
                                digit_size, s_y + 2);

    if (digit_size != last_digit_size) {
      int offset = N_digits * dim;
      int digit_idx = N_digits;
      for (int i = offset; i < num_digits * dim; i++) {
        spec[i] = new uint64_t[last_N];
        digits_exp[i] &= LUT_out_mask;
        for (int j = 0; j < last_N; j++) {
          int idx = (x_digits[i] + j) & last_digit_mask;
          spec[i][j] = (lookup_neg_exp(idx, s_x - digit_size * digit_idx, s_y) -
                        digits_exp[i]) &
                       LUT_out_mask;
        }
      }
      aux->lookup_table<uint64_t>(spec + offset, nullptr, nullptr, dim,
                                  last_digit_size, s_y + 2);
    }

    for (int i = 0; i < num_digits * dim; i++)
      delete[] spec[i];
    delete[] spec;
  } else {
    aux->lookup_table<uint64_t>(nullptr, x_digits, digits_exp, N_digits * dim,
                                digit_size, s_y + 2);
    if (digit_size != last_digit_size) {
      int offset = N_digits * dim;
      aux->lookup_table<uint64_t>(nullptr, x_digits + offset,
                                  digits_exp + offset, dim, last_digit_size,
                                  s_y + 2);
    }
    for (int i = 0; i < num_digits * dim; i++) {
      digits_exp[i] &= LUT_out_mask;
    }
  }

  uint8_t *zero_shares = new uint8_t[dim];
  for (int i = 0; i < dim; i++) {
    zero_shares[i] = 0;
  }
  for (int i = 1; i < num_digits; i *= 2) {
    for (int j = 0; j < num_digits and j + i < num_digits; j += 2 * i) {
      mult->hadamard_product(dim, digits_exp + j * dim,
                             digits_exp + (j + i) * dim, digits_exp + j * dim,
                             s_y + 2, s_y + 2, 2 * s_y + 2, false, false,
                             MultMode::None, zero_shares, zero_shares);
      trunc->truncate_and_reduce(dim, digits_exp + j * dim,
                                 digits_exp + j * dim, s_y, 2 * s_y + 2);
    }
  }
  xt->z_extend(dim, digits_exp, y, s_y + 2, bw_y, zero_shares);

  delete[] x_digits;
  delete[] tmp_1;
  delete[] digits_exp;
  delete[] zero_shares;
}

void MathFunctions::sigmoid(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                            int32_t bw_y, int32_t s_x, int32_t s_y) {
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));
  uint64_t mask_exp_out = ((s_y + 2) == 64 ? -1 : ((1ULL << (s_y + 2)) - 1));
  uint64_t mask_den = ((s_y + 2) == 64 ? -1 : ((1ULL << (s_y + 2)) - 1));
  uint8_t *zero_shares = new uint8_t[dim];
  for (int i = 0; i < dim; i++) {
    zero_shares[i] = 0;
  }

  uint8_t *msb_x = new uint8_t[dim];
  aux->MSB(x, msb_x, dim, bw_x);

  // neg_x = -x + msb_x * (2x) (neg_x is always negative)
  uint64_t *tmp_1 = new uint64_t[dim];
  uint64_t *tmp_2 = new uint64_t[dim];
  for (int i = 0; i < dim; i++) {
    tmp_1[i] = (-1 * x[i]) & mask_x;
    tmp_2[i] = (2 * x[i]) & mask_x;
  }
  uint64_t *neg_x = new uint64_t[dim];
  aux->multiplexer(msb_x, tmp_2, neg_x, dim, bw_x, bw_x);
  for (int i = 0; i < dim; i++) {
    neg_x[i] = (neg_x[i] + tmp_1[i]) & mask_x;
  }

  // den = tmp_1 = 1 + exp_neg_x
  uint64_t *exp_neg_x = new uint64_t[dim];
  lookup_table_exp(dim, neg_x, exp_neg_x, bw_x, s_y + 2, s_x, s_y);
  for (int i = 0; i < dim; i++) {
    tmp_1[i] =
        ((party == ALICE ? (1ULL << s_y) : 0) + exp_neg_x[i]) & mask_exp_out;
  }
  // den can't be 2^{s_y+1}, so 1 is subtracted if msb_den is 1
  uint8_t *msb_den = new uint8_t[dim];
  aux->MSB(tmp_1, msb_den, dim, s_y + 2);
  aux->B2A(msb_den, tmp_2, dim, s_y + 2);
  // den = tmp_1 = den - msb_den
  // tmp_2 = 1 (all 1 vector)
  for (int i = 0; i < dim; i++) {
    tmp_1[i] = (tmp_1[i] - tmp_2[i]) & mask_den;
    tmp_2[i] = (party == ALICE ? 1 : 0);
  }
  uint64_t *sig_neg_x = new uint64_t[dim];
  // sig_neg_x = 1/(1 + exp_neg_x)
  div(dim, tmp_2, tmp_1, sig_neg_x, 2, s_y + 2, s_y + 2, 0, s_y, s_y, true,
      false);

  // tmp_2 = num = 1 + msb_x * (exp_neg_x - 1)
  for (int i = 0; i < dim; i++) {
    tmp_1[i] = (exp_neg_x[i] - (party == ALICE ? 1ULL << s_y : 0)) & mask_den;
  }
  aux->multiplexer(msb_x, tmp_1, tmp_2, dim, s_y + 2, s_y + 2);
  for (int i = 0; i < dim; i++) {
    tmp_2[i] = (tmp_2[i] + (party == ALICE ? 1ULL << s_y : 0)) & mask_den;
  }

  mult->hadamard_product(dim, tmp_2, sig_neg_x, tmp_1, s_y + 2, s_y + 2,
                         2 * s_y + 2, false, false, MultMode::None, zero_shares,
                         zero_shares);
  trunc->truncate_and_reduce(dim, tmp_1, tmp_2, s_y, 2 * s_y + 2);

  if (bw_y <= (s_y + 2)) {
    for (int i = 0; i < dim; i++) {
      y[i] = tmp_2[i] & mask_y;
    }
  } else {
    xt->z_extend(dim, tmp_2, y, s_y + 2, bw_y, zero_shares);
  }

  delete[] msb_x;
  delete[] tmp_1;
  delete[] tmp_2;
  delete[] neg_x;
  delete[] exp_neg_x;
  delete[] msb_den;
  delete[] sig_neg_x;
}

void MathFunctions::tanh(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                         int32_t bw_y, int32_t s_x, int32_t s_y) {
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));
  uint64_t mask_exp_out = ((s_y + 3) == 64 ? -1 : ((1ULL << (s_y + 3)) - 1));
  uint64_t mask_den = ((s_y + 2) == 64 ? -1 : ((1ULL << (s_y + 2)) - 1));

  uint8_t *msb_x = new uint8_t[dim];
  aux->MSB(x, msb_x, dim, bw_x);

  // neg_x = -x + msb_x * (2x) (neg_x is always negative)
  uint64_t *tmp_1 = new uint64_t[dim];
  uint64_t *tmp_2 = new uint64_t[dim];
  for (int i = 0; i < dim; i++) {
    tmp_1[i] = (-1 * x[i]) & mask_x;
    tmp_2[i] = (2 * x[i]) & mask_x;
  }
  uint64_t *neg_x = new uint64_t[dim];
  aux->multiplexer(msb_x, tmp_2, neg_x, dim, bw_x, bw_x);
  for (int i = 0; i < dim; i++) {
    neg_x[i] = (neg_x[i] + tmp_1[i]) & mask_x;
  }

  uint64_t *exp_neg_2x = new uint64_t[dim];
  // scale of neg_x is reduced by 1 instead of mulitplying it with 2 to get
  // neg_2x
  lookup_table_exp(dim, neg_x, exp_neg_2x, bw_x, s_y + 2, s_x - 1, s_y);
  // den = tmp_1 = 1 + exp_neg_2x
  for (int i = 0; i < dim; i++) {
    tmp_1[i] =
        ((party == ALICE ? (1ULL << s_y) : 0) + exp_neg_2x[i]) & mask_exp_out;
  }
  // den can't be 2^{s_y+1}, so 1 is subtracted if msb_den is 1
  uint8_t *msb_den = new uint8_t[dim];
  aux->MSB(tmp_1, msb_den, dim, s_y + 2);
  aux->B2A(msb_den, tmp_2, dim, s_y + 2);
  // den = tmp_1 = den - msb_den
  // num = tmp_2 = 1 - exp_neg_2x
  for (int i = 0; i < dim; i++) {
    tmp_1[i] = (tmp_1[i] - tmp_2[i]) & mask_den;
    tmp_2[i] =
        ((party == ALICE ? (1ULL << s_y) : 0) - exp_neg_2x[i]) & mask_den;
  }
  uint64_t *tanh_neg_x = new uint64_t[dim];
  // tanh_neg_x = (1 - exp_neg_2x)/(1 + exp_neg_2x)
  div(dim, tmp_2, tmp_1, tanh_neg_x, s_y + 2, s_y + 2, s_y + 2, s_y, s_y, s_y,
      true, false);

  // tmp_2 = tanh_neg_x + msb_x * (-2 * tanh_neg_x)
  // tmp_1 = -2 * tanh_neg_x
  for (int i = 0; i < dim; i++) {
    tmp_1[i] = (-2 * tanh_neg_x[i]) & mask_exp_out;
  }
  aux->multiplexer(msb_x, tmp_1, tmp_2, dim, s_y + 2, s_y + 2);
  for (int i = 0; i < dim; i++) {
    tmp_2[i] = (tmp_2[i] + tanh_neg_x[i]) & mask_exp_out;
  }

  if (bw_y <= (s_y + 2)) {
    for (int i = 0; i < dim; i++) {
      y[i] = tmp_2[i] & mask_y;
    }
  } else {
    xt->s_extend(dim, tmp_2, y, s_y + 2, bw_y, msb_x);
  }

  delete[] msb_x;
  delete[] tmp_1;
  delete[] tmp_2;
  delete[] neg_x;
  delete[] exp_neg_2x;
  delete[] msb_den;
  delete[] tanh_neg_x;
}

uint64_t lookup_sqrt(int32_t index, int32_t m, int32_t exp_parity) {
  int32_t k = 1 << m;
  double u = (1.0 + (double(index) / double(k))) * (1 << exp_parity);
  double Y = 1.0 / sqrt(u);
  int32_t scale = m + SQRT_LOOKUP_SCALE;
  uint64_t val = (Y * (1ULL << scale));
  return val;
}

void MathFunctions::sqrt(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                         int32_t bw_y, int32_t s_x, int32_t s_y, bool inverse) {
  int32_t m, iters;
  if (s_y <= 14) {
    m = ceil(s_y / 2.0);
    iters = 1;
  } else {
    m = ceil((ceil(s_y / 2.0)) / 2.0);
    iters = 2;
  }
  assert(m <= KKOT_LIMIT);
  int32_t bw_adjust = bw_x - 1;
  uint64_t mask_adjust = (bw_adjust == 64 ? -1 : ((1ULL << bw_adjust) - 1));
  // MSB is always 0, thus, not including it
  uint8_t *msnzb_vector_bool = new uint8_t[dim * (bw_x - 1)];
  uint64_t *msnzb_vector = new uint64_t[dim * (bw_x - 1)];
  aux->msnzb_one_hot(x, msnzb_vector_bool, bw_x - 1, dim);
  aux->B2A(msnzb_vector_bool, msnzb_vector, dim * (bw_x - 1), bw_x - 1);
  uint64_t *adjust = new uint64_t[dim];
  uint8_t *exp_parity = new uint8_t[dim];
  for (int i = 0; i < dim; i++) {
    adjust[i] = 0;
    exp_parity[i] = 0;
    for (int j = 0; j < (bw_x - 1); j++) {
      adjust[i] += (1ULL << (bw_x - 2 - j)) * msnzb_vector[i * (bw_x - 1) + j];
      if (((j - s_x) & 1)) {
        exp_parity[i] ^= msnzb_vector_bool[i * (bw_x - 1) + j];
      }
    }
    adjust[i] &= mask_adjust;
  }
  // adjusted_x: bw = bw_x - 1, scale = bw_x - 2
  uint64_t *adjusted_x = new uint64_t[dim];
  mult->hadamard_product(dim, x, adjust, adjusted_x, bw_x - 1, bw_x - 1,
                         bw_x - 1, false, false, MultMode::None);
  // Top m bits of adjusted_x excluding MSB, which is always 1
  // adjusted_x_m: bw = m, scale = m
  uint64_t *adjusted_x_m = new uint64_t[dim];
  trunc->truncate_and_reduce(dim, adjusted_x, adjusted_x_m, bw_x - m - 2,
                             bw_x - 2);

  // m + 1 bits will be input to the lookup table
  int32_t M = (1LL << (m + 1));
  uint64_t Y_mask = (1ULL << (m + SQRT_LOOKUP_SCALE + 1)) - 1;
  uint64_t m_mask = (1ULL << m) - 1;
  // Y: bw = m + SQRT_LOOKUP_SCALE + 1, scale = m + SQRT_LOOKUP_SCALE
  uint64_t *Y = new uint64_t[dim];
  if (party == ALICE) {
    uint64_t **spec;
    spec = new uint64_t *[dim];
    PRG128 prg;
    prg.random_data(Y, dim * sizeof(uint64_t));
    for (int i = 0; i < dim; i++) {
      spec[i] = new uint64_t[M];
      Y[i] &= Y_mask;
      for (int j = 0; j < M; j++) {
        // j = exp_parity || (adjusted_x_m) (LSB -> MSB)
        int32_t idx = (adjusted_x_m[i] + (j >> 1)) & m_mask;
        int32_t exp_parity_val = (exp_parity[i] ^ (j & 1));
        spec[i][j] = (lookup_sqrt(idx, m, exp_parity_val) - Y[i]) & Y_mask;
      }
    }
    aux->lookup_table<uint64_t>(spec, nullptr, nullptr, dim, m + 1,
                                m + SQRT_LOOKUP_SCALE + 1);

    for (int i = 0; i < dim; i++)
      delete[] spec[i];
    delete[] spec;
  } else {
    // lut_in = exp_parity || adjusted_x_m
    uint64_t *lut_in = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
      lut_in[i] = ((adjusted_x_m[i] & m_mask) << 1) | (exp_parity[i] & 1);
    }
    aux->lookup_table<uint64_t>(nullptr, lut_in, Y, dim, m + 1,
                                m + SQRT_LOOKUP_SCALE + 1);

    delete[] lut_in;
  }
  // X = (exp_parity ? 2 * adjusted_x : adjusted_x); scale = bw_x - 2
  // X: bw = bw_x
  uint64_t *X = new uint64_t[dim];
  uint64_t *tmp_1 = new uint64_t[dim];
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  xt->z_extend(dim, adjusted_x, X, bw_x - 1, bw_x);
  aux->multiplexer(exp_parity, X, tmp_1, dim, bw_x, bw_x);
  for (int i = 0; i < dim; i++) {
    X[i] = (X[i] + tmp_1[i]) & mask_x;
  }

  uint64_t *x_prev = new uint64_t[dim];
  if (inverse) {
    // x_prev : bw = m + SQRT_LOOKUP_SCALE + 1, scale = m + SQRT_LOOKUP_SCALE
    // x_prev \approx 0.5 < 1/sqrt(X) < 1
    memcpy(x_prev, Y, dim * sizeof(uint64_t));
  } else {
    // x_prev : bw = s_y + 1, scale = s_y
    // x_prev \approx 1 < sqrt(X) < 2
    mult->hadamard_product(dim, X, Y, tmp_1, bw_x, m + SQRT_LOOKUP_SCALE + 1,
                           bw_x - 1 + m + SQRT_LOOKUP_SCALE, false, false,
                           MultMode::None);
    assert((bw_x - 2 + m + SQRT_LOOKUP_SCALE) >= s_y);
    trunc->truncate_and_reduce(dim, tmp_1, x_prev,
                               bw_x - 2 + m + SQRT_LOOKUP_SCALE - s_y,
                               bw_x - 1 + m + SQRT_LOOKUP_SCALE);
  }
  // b_prev: bw = s_y + 2, scale = s_y
  uint64_t *b_prev = new uint64_t[dim];
  assert((bw_x - 2) >= s_y);
  trunc->truncate_and_reduce(dim, X, b_prev, bw_x - 2 - s_y, bw_x);
  // Y_prev: bw = m + SQRT_LOOKUP_SCALE + 1, scale = m + SQRT_LOOKUP_SCALE
  uint64_t *Y_prev = new uint64_t[dim];
  memcpy(Y_prev, Y, dim * sizeof(uint64_t));

  uint64_t b_mask = (1ULL << (s_y + 2)) - 1;
  uint64_t *x_curr = new uint64_t[dim];
  uint64_t *b_curr = new uint64_t[dim];
  uint64_t *Y_curr = new uint64_t[dim];
  uint64_t *Y_sq = new uint64_t[dim];
  for (int i = 0; i < iters; i++) {
    if (i == 0) {
      // Y_sq: bw = 2m + 2SQRT_LOOKUP_SCALE + 1, scale = 2m + 2SQRT_LOOKUP_SCALE
      mult->hadamard_product(
          dim, Y_prev, Y_prev, Y_sq, m + SQRT_LOOKUP_SCALE + 1,
          m + SQRT_LOOKUP_SCALE + 1, 2 * m + 2 * SQRT_LOOKUP_SCALE + 1, false,
          false, MultMode::None);
      // tmp_1: bw = 2m+2SQRT_LOOKUP_SCALE+s_y+2, scale =
      // 2m+2SQRT_LOOKUP_SCALE+s_y
      mult->hadamard_product(dim, Y_sq, b_prev, tmp_1,
                             2 * m + 2 * SQRT_LOOKUP_SCALE + 1, s_y + 2,
                             2 * m + 2 * SQRT_LOOKUP_SCALE + s_y + 2, false,
                             false, MultMode::None);
      // b_curr: bw = s_y + 2, scale = s_y
      trunc->truncate_and_reduce(dim, tmp_1, b_curr,
                                 2 * m + 2 * SQRT_LOOKUP_SCALE,
                                 2 * m + 2 * SQRT_LOOKUP_SCALE + s_y + 2);
    } else {
      // tmp_1: bw = 2*s_y + 3, scale = 2*s_y + 2
      mult->hadamard_product(dim, Y_prev, Y_prev, tmp_1, s_y + 2, s_y + 2,
                             2 * s_y + 3, false, false, MultMode::None);
      // Y_sq: bw = s_y + 1, scale = s_y
      trunc->truncate_and_reduce(dim, tmp_1, Y_sq, s_y + 2, 2 * s_y + 3);
      // tmp_1: bw = 2s_y + 2, scale = 2s_y
      mult->hadamard_product(dim, Y_sq, b_prev, tmp_1, s_y + 1, s_y + 2,
                             2 * s_y + 2, false, false, MultMode::None);
      // b_curr: bw = s_y + 2, scale = s_y
      trunc->truncate_and_reduce(dim, tmp_1, b_curr, s_y, 2 * s_y + 2);
    }
    for (int j = 0; j < dim; j++) {
      // Y_curr: bw = s_y + 2, scale = s_y + 1
      // Y_curr = (3 - b_curr)/2
      Y_curr[j] = ((party == ALICE ? (3ULL << s_y) : 0) - b_curr[j]) & b_mask;
    }
    if (inverse && (i == 0)) {
      // tmp_1: bw = s_y+m+SQRT_LOOKUP_SCALE+2, scale =
      // s_y+m+SQRT_LOOKUP_SCALE+1
      mult->hadamard_product(
          dim, x_prev, Y_curr, tmp_1, m + SQRT_LOOKUP_SCALE + 1, s_y + 2,
          s_y + m + SQRT_LOOKUP_SCALE + 2, false, false, MultMode::None);
      // x_curr: bw = s_y + 1, scale = s_y
      trunc->truncate_and_reduce(dim, tmp_1, x_curr, m + SQRT_LOOKUP_SCALE + 1,
                                 s_y + m + SQRT_LOOKUP_SCALE + 2);
    } else {
      // tmp_1: bw = 2*s_y + 2, scale = 2s_y + 1
      mult->hadamard_product(dim, x_prev, Y_curr, tmp_1, s_y + 1, s_y + 2,
                             2 * s_y + 2, false, false, MultMode::None);
      // x_curr: bw = s_y + 1, scale = s_y
      trunc->truncate_and_reduce(dim, tmp_1, x_curr, s_y + 1, 2 * s_y + 2);
    }
    memcpy(x_prev, x_curr, dim * sizeof(uint64_t));
    memcpy(b_prev, b_curr, dim * sizeof(uint64_t));
    memcpy(Y_prev, Y_curr, dim * sizeof(uint64_t));
  }

  int32_t bw_sqrt_adjust = bw_x / 2;
  uint64_t mask_sqrt_adjust =
      (bw_sqrt_adjust == 64 ? -1 : ((1ULL << bw_sqrt_adjust) - 1));
  uint64_t *sqrt_adjust = new uint64_t[dim];
  int32_t sqrt_adjust_scale =
      (inverse ? floor((bw_x - 1 - s_x) / 2.0) : floor((s_x + 1) / 2.0));
  for (int i = 0; i < dim; i++) {
    sqrt_adjust[i] = 0;
    for (int j = 0; j < (bw_x - 1); j++) {
      if (inverse) {
        sqrt_adjust[i] +=
            (1ULL << int(floor((s_x - j + 1) / 2.0) + sqrt_adjust_scale)) *
            msnzb_vector[i * (bw_x - 1) + j];
      } else {
        sqrt_adjust[i] +=
            (1ULL << int(floor((j - s_x) / 2.0) + sqrt_adjust_scale)) *
            msnzb_vector[i * (bw_x - 1) + j];
      }
    }
    sqrt_adjust[i] &= mask_sqrt_adjust;
  }
  if (iters > 0 || (!inverse)) {
    // tmp_1: bw = s_y + 1 + bw_sqrt_adjust, scale = s_y + sqrt_adjust_scale
    mult->hadamard_product(dim, x_prev, sqrt_adjust, tmp_1, s_y + 1,
                           bw_sqrt_adjust, s_y + 1 + bw_sqrt_adjust, false,
                           false, MultMode::None);
    // x_curr: bw = s_y + floor(bw_x/2) + 1 - ceil(s_x/2), scale = s_y
    trunc->truncate_and_reduce(dim, tmp_1, x_prev, sqrt_adjust_scale,
                               s_y + 1 + bw_sqrt_adjust);
    if (bw_y > (s_y + 1 + bw_sqrt_adjust - sqrt_adjust_scale)) {
      xt->z_extend(dim, x_prev, y, s_y + 1 + bw_sqrt_adjust - sqrt_adjust_scale,
                   bw_y);
    } else {
      aux->reduce(dim, x_prev, y, s_y + 1 + bw_sqrt_adjust - sqrt_adjust_scale,
                  bw_y);
    }
  } else {
    // tmp_1: bw = m + SQRT_LOOKUP_SCALE + 1 + bw_sqrt_adjust,
    //        scale = m + SQRT_LOOKUP_SCALE + sqrt_adjust_scale
    mult->hadamard_product(dim, x_prev, sqrt_adjust, tmp_1,
                           m + SQRT_LOOKUP_SCALE + 1, bw_sqrt_adjust,
                           m + SQRT_LOOKUP_SCALE + 1 + bw_sqrt_adjust, false,
                           false, MultMode::None);
    // x_curr: bw = m + floor(bw_x/2) + 1 - ceil(s_x/2), scale = m
    // If iters == 0, we know s_y = m
    trunc->truncate_and_reduce(dim, tmp_1, x_prev,
                               sqrt_adjust_scale + SQRT_LOOKUP_SCALE,
                               m + SQRT_LOOKUP_SCALE + 1 + bw_sqrt_adjust);
    if (bw_y > (m + 1 + bw_sqrt_adjust - sqrt_adjust_scale)) {
      xt->z_extend(dim, x_prev, y, m + 1 + bw_sqrt_adjust - sqrt_adjust_scale,
                   bw_y);
    } else {
      aux->reduce(dim, x_prev, y, m + 1 + bw_sqrt_adjust - sqrt_adjust_scale,
                  bw_y);
    }
  }

  delete[] msnzb_vector_bool;
  delete[] msnzb_vector;
  delete[] adjust;
  delete[] exp_parity;
  delete[] adjusted_x;
  delete[] X;
  delete[] tmp_1;
  delete[] x_prev;
  delete[] b_prev;
  delete[] Y_prev;
  delete[] x_curr;
  delete[] b_curr;
  delete[] Y_curr;
  delete[] Y_sq;
  delete[] sqrt_adjust;

  return;
}

void MathFunctions::ReLU(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                         uint64_t six) {
  bool six_comparison = false;
  if (six != 0)
    six_comparison = true;
  int32_t size = (six_comparison ? 2 * dim : dim);

  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

  uint64_t *tmp = new uint64_t[size];
  uint8_t *tmp_msb = new uint8_t[size];
  memcpy(tmp, x, dim * sizeof(uint64_t));
  if (six_comparison) {
    for (int i = 0; i < dim; i++) {
      tmp[dim + i] = (party == ALICE ? x[i] - six : x[i]) & mask_x;
    }
  }
  aux->MSB(tmp, tmp_msb, size, bw_x);
  for (int i = 0; i < size; i++) {
    if (party == ALICE) {
      tmp_msb[i] = tmp_msb[i] ^ 1;
    }
  }
  if (six_comparison)
    aux->AND(tmp_msb, tmp_msb + dim, tmp_msb + dim, dim);

  aux->multiplexer(tmp_msb, tmp, tmp, size, bw_x, bw_x);

  memcpy(y, tmp, dim * sizeof(uint64_t));
  if (six_comparison) {
    for (int i = 0; i < dim; i++) {
      y[i] = (y[i] - tmp[i + dim]) & mask_x;
    }
  }

  delete[] tmp;
  delete[] tmp_msb;
}

/*
  cryptoLp
*/

void MathFunctions::ABS(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x) {
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint8_t *msb = new uint8_t[dim];

  aux->MSB(x, msb, dim, bw_x);

  for (int i = 0; i < dim; i++) {
    msb[i] = msb[i] & 1;
    x[i] = x[i] & mask_x;
  }

  aux->ABS_multiplexer(msb, x, y, dim, bw_x, bw_x);

  delete[] msb;
}

void MathFunctions::diff(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x) {
  assert(bw_x <= 64);
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t *s12 = new uint64_t[dim]; // share of 'x_b' from party 'b' to '1-b'
  uint64_t *s21 = new uint64_t[dim]; // share of 'x_{1-b}' from party '1-b' to 'b'
  uint64_t *r = new uint64_t[dim]; // random mask

  PRG128 prg;
  prg.random_data(r, dim * sizeof(uint64_t));
  for (int i = 0; i < dim; i++) {
    r[i] &= mask_x;
    s12[i] = (party == ALICE ? r[i] - x[i] : x[i] + r[i]) & mask_x;
  }
  iopack->io->send_data(s12, dim * sizeof(uint64_t));
  iopack->io->recv_data(s21, dim * sizeof(uint64_t));

  for (int i = 0; i < dim; i++) {
    y[i] = (s21[i] - r[i]) & mask_x;
  }

  delete[] r;
  delete[] s12;
  delete[] s21;
}

void MathFunctions::millionaire(uint64_t *x, uint8_t *y, int32_t size,
                               int32_t bw_x) {
  assert(bw_x <= 64); // l <= 64
  uint64_t mask = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

  uint64_t *tmp_x = new uint64_t[size];
  for (int i = 0; i < size; i++) {
    tmp_x[i] = x[i] & mask;
  }
  aux->mill->compare(y, tmp_x, size, bw_x, false); // computing less_than

  delete[] tmp_x;
}

void MathFunctions::millionaire_ext_naive(uint64_t *x, uint8_t *y, int32_t size, 
                                          int32_t bw_x) {
  assert(bw_x <= 64); // l <= 64
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t pow_x = (bw_x == 64 ? 0ULL : (1ULL << bw_x));
  uint64_t half_pow_x = pow_x >> 1;

  uint64_t *d = new uint64_t[size];
  uint8_t *msb_d = new uint8_t[size];
  uint8_t *lt_d = new uint8_t[size];

  diff(size, x, d, bw_x);

  aux->MSB(d, msb_d, size, bw_x);

  for (int i = 0; i < size; i++) {
    lt_d[i] = (party == ALICE ? msb_d[i] ^ 1 : msb_d[i]);
  }

  uint8_t *lt = new uint8_t[size];
  
  for (int i = 0; i < size; i++) {
    lt[i] = (x[i] < half_pow_x);
  }

  uint8_t *w_xnot = new uint8_t[size];
  uint8_t *xnot = new uint8_t[size];
  uint8_t *w = new uint8_t[size];
  for (int i = 0; i < size; i++) {
    xnot[i] = (party == ALICE ? lt[i] ^ 1 : 0);
    w[i] = (party == ALICE ? 0 : lt[i]);
  }
  aux->AND(xnot, w, w_xnot, size);

  uint8_t *wx = new uint8_t[size];
  uint8_t *temp = new uint8_t[size];
  for (int i = 0; i < size; i++) {
    temp[i] = (party == ALICE ? lt[i] : 0);
  }
  aux->AND(w, temp, wx, size);

  uint8_t *wnot_xnot_ynot = new uint8_t[size];
  uint8_t *wnot_xnot = new uint8_t[size];
  uint8_t *wnot = new uint8_t[size];
  uint8_t *ynot = new uint8_t[size];
  for (int i = 0; i < size; i++) {
    wnot[i] = (party == ALICE ? 0 : lt[i] ^ 1);
    ynot[i] = (party == ALICE ? lt_d[i] ^ 1 : lt_d[i]);
  }
  aux->AND(wnot, xnot, wnot_xnot, size);
  aux->AND(wnot_xnot, ynot, wnot_xnot_ynot, size);

  uint8_t *wx_ynot = new uint8_t[size];
  for (int i = 0; i < size; i++) {
    ynot[i] = (party == ALICE ? lt_d[i] ^ 1 :lt_d[i]);
  }
  aux->AND(wx, ynot, wx_ynot, size);

  for (int i = 0; i < size; i++) {
    y[i] = (w_xnot[i] ^ wnot_xnot_ynot[i] ^ wx_ynot[i]);
  }

  delete[] d;
  delete[] msb_d;
  delete[] lt_d;
  delete[] lt;
  delete[] w_xnot;
  delete[] xnot;
  delete[] w;
  delete[] wx;
  delete[] temp;
  delete[] wnot_xnot_ynot;
  delete[] wnot_xnot;
  delete[] wnot;
  delete[] ynot;
  delete[] wx_ynot;
}

void MathFunctions::millionaire_ext_optimized(uint64_t *x, uint8_t *y, int32_t size, 
                                              int32_t bw_x) {
  assert(bw_x <= 64); // l <= 64
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t pow_x = (bw_x == 64 ? 0ULL : (1ULL << bw_x));
  uint64_t half_pow_x = pow_x >> 1;

  uint64_t *d = new uint64_t[size];
  uint8_t *msb_d = new uint8_t[size];
  uint8_t *lt = new uint8_t[size];
  uint8_t *lt_d = new uint8_t[size];
  uint8_t *xy = new uint8_t[size];
  uint8_t *temp_x = new uint8_t[size];
  uint8_t *w = new uint8_t[size];

  diff(size, x, d, bw_x);

  aux->MSB(d, msb_d, size, bw_x);
  for (int i = 0; i < size; i++) {
    lt_d[i] = (party == ALICE ? msb_d[i] ^ 1 : msb_d[i]);
    lt[i] = (x[i] < half_pow_x);
    temp_x[i] = (party == ALICE ? lt[i] : 0);
    w[i] = (party == ALICE ? 0 : lt[i]);
  }
  aux->AND(temp_x, lt_d, xy, size);

  // temp1 = w(x + y - 2xy)
  uint8_t *temp1 = new uint8_t[size];
  uint8_t *temp2 = new uint8_t[size];
  uint8_t *temp3 = new uint8_t[size];
  for (int i = 0; i < size; i++) {
    temp1[i] = (temp_x[i] + lt_d[i] - 2 * xy[i]) & 1;
  }
  aux->AND(w, temp1, temp2, size);

  for (int i = 0; i < size; i++) {
    temp3[i] = (temp2[i] - lt_d[i] - temp_x[i] + xy[i]) & 1;
    y[i] = (party == ALICE ? temp3[i] : (temp3[i] + 1) & 1);
  }
}

void MathFunctions::millionaire_ext(uint64_t *x, uint8_t *y, int32_t size, 
                                    int32_t bw_x) {
  assert(bw_x <= 64); // l <= 64
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

  int beta = 2; // m
  int beta_pow = 1 << beta;
  uint8_t mask_beta = (1 << beta) - 1;

  TripleGenerator *triple_gen;
  triple_gen = new TripleGenerator(party, iopack, otpack);

  int num_digits = ceil((double)bw_x / beta); // q
  uint8_t *digits = new uint8_t[num_digits * size];
  uint8_t *lt = new uint8_t[num_digits * size];
  uint8_t *eq = new uint8_t[num_digits * size];
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < num_digits; j++) {
      digits[i * num_digits + j] = (uint8_t)(x[i] >> j * beta) & mask_beta;
    }
  }
  
  if (party == sci::ALICE) {
    uint8_t * *ot_messages;
    ot_messages = new uint8_t *[num_digits * size];
    for (int i = 0; i < num_digits * size; i++)
      ot_messages[i] = new uint8_t[beta_pow];
    triple_gen->prg->random_bool((bool *)lt, num_digits * size);
    triple_gen->prg->random_bool((bool *)eq, num_digits * size);
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < num_digits; j++) {
        aux->mill->set_leaf_ot_messages(ot_messages[i * num_digits + j],
                                        digits[i * num_digits + j], beta_pow,
                                        lt[i * num_digits + j], eq[i * num_digits + j],
                                        false);
      }
    }
    otpack->kkot[beta - 1]->send(ot_messages, size * num_digits, 2);
    // Cleanup
    for (int i = 0; i < num_digits * size; i++)
      delete[] ot_messages[i];
    delete[] ot_messages;
  } else {
    otpack->kkot[beta - 1]->recv(lt, digits, size * num_digits, 2);
    for (int i = 0; i < num_digits * size; i++) {
      eq[i] = lt[i] & 1;
      lt[i] >>= 1;
    }
  }

  uint8_t *lt_AND = new uint8_t[(num_digits - 1) * size];
  uint8_t *eq_AND = new uint8_t[(num_digits - 1) * size];
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < num_digits; j++) {
      if (j != num_digits - 1) {
        lt_AND[i * (num_digits - 1) + j] = lt[i * num_digits + j];
        eq_AND[i * (num_digits - 1) + j] = eq[i * num_digits + (j + 1)];
      }
    }
  }

  uint8_t *curr =  new uint8_t[1];
  uint8_t *pre =  new uint8_t[1];
  uint8_t * *eq_cum;
  uint8_t *eq_cum_cp = new uint8_t[(num_digits - 1) * size];
  eq_cum = new uint8_t *[(num_digits - 1) * size];
  for (int i = 0; i < size; i++) {
    for (int j = num_digits - 2; j >= 0; j--) {
      eq_cum[i * (num_digits - 1) + j] = new uint8_t[1];
      if (j != num_digits - 2) {
        curr[0] = eq_AND[i * (num_digits - 1) + j];
        pre[0] = eq_cum[i * (num_digits - 1) + (j + 1)][0];
        aux->AND(curr, pre, eq_cum[i * (num_digits - 1) + j], 1);
        eq_cum_cp[i * (num_digits - 1) + j] = eq_cum[i * (num_digits - 1) + j][0];
      } else {
        eq_cum[i * (num_digits - 1) + j][0] = eq_AND[i * (num_digits - 1) + j];
        eq_cum_cp[i * (num_digits - 1) + j] = eq_cum[i * (num_digits - 1) + j][0];
      }
    }
  }

  uint8_t *lt_eq = new uint8_t[(num_digits - 1) * size];
  aux->AND(lt_AND, eq_cum_cp, lt_eq, (num_digits - 1) * size);

  // // Test And
  // for (int i = 0; i < size; i++) {
  //   for (int j = 0; j < num_digits - 1; j++) {
  //     y[i * (num_digits - 1) + j] = lt_eq[i * (num_digits - 1) + j];
  //   }
  // }

  uint8_t *res = new uint8_t[size](); // Initialized to 0
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < num_digits; j++) {
      if (j == num_digits - 1) {
        res[i] ^= lt[i * num_digits + j];
      } else {
        res[i] ^= lt_eq[i * (num_digits - 1) + j];
      }
    }
  }

  for (int i = 0; i < size; i++) {
    y[i] = res[i];
  }

  // Cleanup
  delete[] digits;
  delete[] lt;
  delete[] eq;
  delete[] res;
  delete[] lt_eq;
  delete triple_gen;
}

void MathFunctions::millionaire_ext_plus(uint64_t *x, uint8_t *y, int32_t size, 
                                         int32_t bw_x) {
  assert(bw_x <= 64); // l <= 64
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

  int beta = 4; // m
  int beta_pow = 1 << beta;
  uint8_t mask_beta = (1 << beta) - 1;

  TripleGenerator *triple_gen;
  triple_gen = new TripleGenerator(party, iopack, otpack);

  int num_digits = ceil((double)bw_x / beta); // q
  uint8_t *digits = new uint8_t[num_digits * size];
  uint8_t *lt = new uint8_t[num_digits * size];
  uint8_t *eq = new uint8_t[num_digits * size];
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < num_digits; j++) {
      digits[i * num_digits + j] = (uint8_t)(x[i] >> j * beta) & mask_beta;
    }
  }
  
  if (party == sci::ALICE) {
    uint8_t * *ot_messages;
    ot_messages = new uint8_t *[num_digits * size];
    for (int i = 0; i < num_digits * size; i++)
      ot_messages[i] = new uint8_t[beta_pow];
    triple_gen->prg->random_bool((bool *)lt, num_digits * size);
    triple_gen->prg->random_bool((bool *)eq, num_digits * size);
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < num_digits; j++) {
        aux->mill->set_leaf_ot_messages(ot_messages[i * num_digits + j],
                                        digits[i * num_digits + j], beta_pow,
                                        lt[i * num_digits + j], eq[i * num_digits + j],
                                        false);
      }
    }
    otpack->kkot[beta - 1]->send(ot_messages, size * num_digits, 2);
    // Cleanup
    for (int i = 0; i < num_digits * size; i++)
      delete[] ot_messages[i];
    delete[] ot_messages;
  } else {
    otpack->kkot[beta - 1]->recv(lt, digits, size * num_digits, 2);
    for (int i = 0; i < num_digits * size; i++) {
      eq[i] = lt[i] & 1;
      lt[i] >>= 1;
    }
  }

  uint8_t * *res;
  res = new uint8_t *[size];
  uint8_t *left = new uint8_t[1];
  uint8_t *right = new uint8_t[1];
  for (int i = 0; i < size; i++) {
    res[i] = new uint8_t[1]();
    for (int j = 0; j < num_digits - 1; j++) {
      right[0] = lt[i * num_digits + j] ^ res[i][0];
      left[0] = eq[i * num_digits + (j + 1)];
      aux->AND(left, right, res[i], 1);
    }
    res[i][0] ^= lt[i * num_digits + (num_digits - 1)];
  }

  for (int i = 0; i < size; i++) {
    y[i] = res[i][0];
  }

  // Cleanup
  delete[] digits;
  delete[] lt;
  delete[] eq;
  for (int i = 0; i < size; i++)
    delete[] res[i];
  delete triple_gen;
}

void MathFunctions::MD(int32_t dim, uint64_t *x, uint64_t *y, uint64_t *z, 
                       int32_t bw_x) {
  assert(bw_x <= 64);
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t *d = new uint64_t[dim];

  for (int i = 0; i < dim; i++) {
    d[i] = (x[i] - y[i]) & mask_x;
  }

  ABS(dim, d, z, bw_x);

  delete[] d;
}

// each party has plaintexts x1 and x2, and they want to compute the dot product
void MathFunctions::elemwise_prod(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x, 
                                  int32_t bw_y, MultMode mode) {
  assert(mode == MultMode::None);
  assert(bw_x <= 64 && bw_y <= 64);
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));

  uint64_t *temp1 = new uint64_t[dim];
  uint64_t *temp2 = new uint64_t[dim];

  for (int i = 0; i < dim; i++) {
    temp1[i] = (party == ALICE) ? x[i] : 0;
    temp2[i] = (party == ALICE) ? 0 : x[i];
  }

  mult->hadamard_product(dim, temp1, temp2, y, bw_x, bw_x, bw_y, false, false, mode);

  delete[] temp1;
  delete[] temp2;
}

void MathFunctions::ED(int32_t dim, uint64_t *x, uint64_t *y, uint64_t *z,
                       int32_t bw_x, int32_t bw_z) {
  assert(bw_x <= 64 && bw_z <= 64);
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t mask_z = (bw_z == 64 ? -1 : ((1ULL << bw_z) - 1));
  uint64_t *d_sqr = new uint64_t[dim];

  uint64_t *d = new uint64_t[dim];

  for (int i = 0; i < dim; i++) {
    d[i] = (x[i] - y[i]) & mask_x;
  }

  // TODO: Can you do sqr faster than mult?
  mult->hadamard_product(dim, d, d, d_sqr, bw_x, bw_x, bw_z, false, false, MultMode::None);

  for (int i = 0; i < dim; i++) {
    z[i] = d_sqr[i];
  }

  // sqrt(dim, d_sqr, z, bw_x, bw_z, 12, 12, false);

  delete[] d_sqr;
}

void MathFunctions::ED_alt(int32_t dim, uint64_t *x, uint64_t *y, uint64_t *z,
                           int32_t bw_x, int32_t bw_z) {
  assert(bw_x <= 64 && bw_z <= 64);
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t mask_z = (bw_z == 64 ? -1 : ((1ULL << bw_z) - 1));
  uint64_t *d_sqr = new uint64_t[dim];
  uint64_t *prod = new uint64_t[dim];

  uint64_t *d = new uint64_t[dim];

  for (int i = 0; i < dim; i++) {
    d[i] = (x[i] - y[i]) & mask_x;
    d_sqr[i] = (d[i] * d[i]) & mask_z;
  }

  elemwise_prod(dim, d, prod, bw_x, bw_z, MultMode::None);

  for (int i = 0; i < dim; i++) {
    z[i] = (d_sqr[i] + 2 * prod[i]) & mask_z;
  }

  // sqrt(dim, d_sqr, y, bw_x, bw_z, 12, 12);

  delete[] d_sqr;
  delete[] prod;
}

// naiive implementation
void MathFunctions::max_naive(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x) {
  assert(bw_x <= 64);
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

  uint64_t *max = new uint64_t[1];
  max[0] = x[0] & mask_x;
  uint64_t *d = new uint64_t[1];
  for (int i = 1; i < dim; i++) {
    d[0] = (max[0] - x[i]) & mask_x;

    uint8_t *msb = new uint8_t[1];
    aux->MSB(d, msb, 1, bw_x);
    msb[0] = party == ALICE ? msb[0] ^ 1 : msb[0];
    aux->multiplexer(msb, d, max, 1, bw_x, bw_x);
    max[0] = (max[0] + x[i]) & mask_x;
  }
  
  y[0] = max[0];

  delete[] d;
  delete[] max;
}

void MathFunctions::max(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x) {
  assert(bw_x <= 64);
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t *d = new uint64_t[dim];

  int n = dim;
  while(n > 1) {
    int new_dim = (n + 1) / 2; // number of nodes in the next level
    uint64_t *d = new uint64_t[new_dim];
    for (int i = 0; i < new_dim; ++i) {
      if (!(n % 2)) { // check if there is any unmatched node
        d[i] = x[2 * i] - x[2 * i + 1];
      } else { // move up the unmatched node
        if (i != new_dim - 1) {
          d[i] = x[2 * i] - x[2 * i + 1];
        } else {
          d[i] = x[2 * i];
        }
      }
    }
    uint8_t *msb = new uint8_t[new_dim];
    aux->MSB(d, msb, new_dim, bw_x);
    for (int i = 0; i < new_dim; ++i) {
      if (party == ALICE) {
        msb[i] = msb[i] ^ 1;
      }
    }
    uint64_t *m = new uint64_t[new_dim];
    aux->multiplexer(msb, d, m, new_dim, bw_x, bw_x);
    for(int i = 0; i < new_dim; ++i) {
      if (!(n % 2)) {
        m[i] = (m[i] + x[2 * i + 1]) & mask_x;
      } else {
        if (i != new_dim - 1) {
          m[i] = (m[i] + x[2 * i + 1]) & mask_x;
        }
      }
      x[i] = m[i];
    }

    n = new_dim;

    if (n == 1) {
      y[0] = m[0];
    }

    delete[] d;
    delete[] msb;
    delete[] m;
  }
}

void MathFunctions::CD_naive_alt(int32_t dim, uint64_t *x, uint64_t *y, uint64_t *z,
                                 int32_t bw_x) {
  ReLURingProtocol<uint64_t> *relu_oracle = new ReLURingProtocol<uint64_t>(party, RING, iopack, 
                                                                           bw_x, 4, otpack);

  assert(bw_x <= 64);
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t *d = new uint64_t[dim];

  for (int i = 0; i < dim; i++) {
    d[i] = (x[i] - y[i]) & mask_x;
  }

  uint64_t *d_abs = new uint64_t[dim];
  ABS(dim, d, d_abs, bw_x);

  uint64_t *max = new uint64_t[1];
  max[0] = d_abs[0] & mask_x;
  uint64_t *w = new uint64_t[1];
  for (int i = 1; i < dim; i++) {
    w[0] = (max[0] - d_abs[i]) & mask_x;
    relu_oracle->relu(max, w, 1);
    max[0] = (max[0] + d_abs[i]) & mask_x;
  }
  
  z[0] = max[0];

  delete[] d;
  delete[] d_abs;
  delete[] max;
  delete[] w;
}

void MathFunctions::CD(int32_t dim, uint64_t *x, uint64_t *y, uint64_t *z,
                       int32_t bw_x) {
  if (dim == 0) return;

  assert(bw_x <= 64);
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t *d = new uint64_t[dim];

  for (int i = 0; i < dim; i++) {
    d[i] = (x[i] - y[i]) & mask_x;
  }

  uint64_t *d_abs = new uint64_t[dim];
  ABS(dim, d, d_abs, bw_x);

  max(dim, d_abs, z, bw_x);

  delete[] d;
  delete[] d_abs;
}

void MathFunctions::adder(int32_t dim, uint64_t *x, uint64_t *y, uint64_t *z,
                          int32_t bw_x) {
  assert(bw_x <= 64);
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t *temp_z = new uint64_t[dim];

  MD(dim, x, y, temp_z, bw_x);

  for (int i = 0; i < dim; i++) {
    z[i] = (-temp_z[i]) & mask_x;
  }

  delete[] temp_z;
}

void MathFunctions::conv(int32_t dim, uint64_t *x, uint64_t *y, uint64_t *z, int32_t bw_x,
                         int32_t bw_y, int32_t bw_z) {
  assert(bw_x <= 64 && bw_y <= 64 && bw_z <= 64);
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
  uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));
  uint64_t mask_z = (bw_z == 64 ? -1 : ((1ULL << bw_z) - 1));

  mult->hadamard_product(dim, x, y, z, bw_x, bw_y, bw_z, false, false, MultMode::None);
}

void MathFunctions::argmin(int32_t dim, uint64_t *x, uint64_t *y, uint64_t *yi,
                           int32_t bw_x) {
  assert(bw_x <= 64);
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

  int n = dim;
  while (n > 1) {
    uint64_t *index = new uint64_t[n](); // initialize to 0
    for (int i = 0; i < n; i++) {
      if (party == ALICE) {
        index[i] = i;
      }
    }
    int new_dim = (n + 1) / 2; // number of nodes in the next level
    uint64_t *w = new uint64_t[new_dim];
    uint64_t *wi = new uint64_t[new_dim];
    for (int i = 0; i < new_dim; ++i) {
      if (!(n % 2)) { // check if there is any unmatched node
        w[i] = (x[2 * i] - x[2 * i + 1]) & mask_x;
        wi[i] = (index[2 * i] - index[2 * i + 1]) & mask_x;
      } else { // move up the unmatched node
        if (i != new_dim - 1) {
          w[i] = (x[2 * i] - x[2 * i + 1]) & mask_x;
          wi[i] = (index[2 * i] - index[2 * i + 1]) & mask_x;
        } else {
          w[i] = x[2 * i];
          wi[i] = index[2 * i];
        }
      }
    }
    uint8_t *msb = new uint8_t[new_dim];
    aux->MSB(w, msb, new_dim, bw_x);
    for (int i = 0; i < new_dim; ++i) {
      if (party == ALICE) {
        msb[i] = msb[i] ^ 1;
      }
    }
    uint64_t *min = new uint64_t[new_dim];
    uint64_t *argmin = new uint64_t[new_dim];
    aux->multiplexer(msb, w, min, new_dim, bw_x, bw_x);
    aux->multiplexer(msb, wi, argmin, new_dim, bw_x, bw_x);
    for(int i = 0; i < new_dim; ++i) {
      if (!(n % 2)) {
        min[i] = (x[2 * i] - min[i]) & mask_x;
        argmin[i] = (index[2 * i] - argmin[i]) & mask_x;
      } else {
        if (i != new_dim - 1) {
          min[i] = (x[2 * i] - min[i]) & mask_x;
          argmin[i] = (index[2 * i] - argmin[i]) & mask_x;
        }
      }
      x[i] = min[i];
    }

    n = new_dim;

    if (n == 1) {
      y[0] = min[0];
      yi[0] = argmin[0];
    }
    delete[] w;
    delete[] msb;
    delete[] min;
  }
};

void MathFunctions::compare_and_swap_naive(uint64_t *x, int32_t i1, int32_t i2, int32_t dim, 
                                           int32_t bw_x) {
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

  uint64_t *w = new uint64_t[1];
  w[0] = (x[i1] - x[i2]) & mask_x;

  uint8_t *msb = new uint8_t[1];
  uint64_t *z = new uint64_t[1];
  aux->MSB(w, msb, 1, bw_x);
  msb[0] = party == ALICE ? msb[0] ^ 1 : msb[0];
  aux->multiplexer(msb, w, z, 1, bw_x, bw_x);
  uint64_t min;
  uint64_t max;
  min = (x[i1] - z[0]) & mask_x;
  max = (x[i2] + z[0]) & mask_x;

  x[i1] = min;
  x[i2] = max;
}

void MathFunctions::compare_and_swap_bitonic(uint64_t *x, int *i1, int *i2, int32_t batch_size, 
                                             int32_t bw_x, bool ascending) {
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

  uint64_t *w = new uint64_t[batch_size];
  for (int i = 0; i < batch_size; i++) {
    w[i] = (x[i1[i]] - x[i2[i]]) & mask_x;
  }

  uint64_t *z = new uint64_t[batch_size];
  uint8_t *msb = new uint8_t[batch_size];
  aux->MSB(w, msb, batch_size, bw_x);
  for (int i = 0; i < batch_size; i++) {
    msb[i] = party == ALICE ? msb[i] ^ 1 : msb[i];
  }
  aux->multiplexer(msb, w, z, batch_size, bw_x, bw_x);

  for (int i = 0; i < batch_size; i++) {
    uint64_t min;
    uint64_t max;
    min = (x[i1[i]] - z[i]) & mask_x;
    max = (x[i2[i]] + z[i]) & mask_x;
    if (ascending) {
      x[i1[i]] = min;
      x[i2[i]] = max;
    } else {
      x[i1[i]] = max;
      x[i2[i]] = min;
    }
  }
}

void MathFunctions::sort_naive(uint64_t *x, int32_t dim, int32_t bw_x) {
  for (int i = 0; i < dim - 1; i++) {
    for (int j = i + 1; j < dim; j++) {
      compare_and_swap_naive(x, i, j, dim, bw_x);
    }
  }
}

void MathFunctions::sort_bitonic(uint64_t *x, int32_t dim, int32_t bw_x) {
  int n = dim;
  for (int k = 2; k <= n; k *= 2) {
    for (int j = k >> 1; j > 0; j /= 2) {
      int *local_indices_i_asc = new int[n/2];
      int *local_indices_ij_asc = new int[n/2];
      int *local_indices_i_desc = new int[n/2];
      int *local_indices_ij_desc = new int[n/2];
      int local_batch_count_asc = 0;
      int local_batch_count_desc = 0;
      for (int i = 0; i < n; i++) {
        int ij = i ^ j;
        if ((ij) > i) {
          if ((i & k) == 0) {
            local_indices_i_asc[local_batch_count_asc] = i;
            local_indices_ij_asc[local_batch_count_asc] = ij;
            local_batch_count_asc++;
          } else {
            local_indices_i_desc[local_batch_count_desc] = i;
            local_indices_ij_desc[local_batch_count_desc] = ij;
            local_batch_count_desc++;
          }
        }
      }
      if (local_batch_count_asc > 0) {
        compare_and_swap_bitonic(x, local_indices_i_asc, local_indices_ij_asc, local_batch_count_asc, bw_x, true);
      }
      if (local_batch_count_desc > 0) {
        compare_and_swap_bitonic(x, local_indices_i_desc, local_indices_ij_desc, local_batch_count_desc, bw_x, false);
      }
    }
  }
}

void MathFunctions::median(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x) {
  assert(bw_x <= 64);
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

  sort_naive(x, dim, bw_x);
  
  if (dim & 1) { // odd
    y[0] = x[(dim - 1) / 2];
  } else { // even
    uint64_t *temp0 = new uint64_t[1];
    uint64_t *temp1 = new uint64_t[1];
    uint64_t *y0 = new uint64_t[1];
    uint64_t *y1 = new uint64_t[1];
    temp0[0] = x[dim / 2 - 1];
    temp1[0] = x[dim / 2];
    trunc->div_pow2(1, temp0, y0, 1, bw_x, true);
    trunc->div_pow2(1, temp1, y1, 1, bw_x, true);
    y[0] = (y0[0] + y1[0]) & mask_x;
    delete[] temp0;
    delete[] temp1;
    delete[] y0;
    delete[] y1;
  }
}

void MathFunctions::K_medians(int32_t max_iter, int32_t dim, uint64_t *x, uint64_t *y, 
                              int32_t bw_x, int32_t K) {
  assert(bw_x <= 64);
  uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

  uint64_t *medians = new uint64_t[K];
  std::vector<uint64_t>* clusters = new std::vector<uint64_t>[K];
  
  // choose K initial medians
  for (int i = 0; i < K; i++) {
    medians[i] = x[i];
  }

  for (int l = 0; l < max_iter; l++) {
    uint64_t *temp = new uint64_t[K];
    uint64_t *dist = new uint64_t[K];
    uint64_t *min = new uint64_t[1];
    uint64_t *min_index = new uint64_t[1];
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < K; j++) {
        temp[j] = x[i];
      }
      MD(K, temp, medians, dist, bw_x);
      argmin(K, dist, min, min_index, bw_x);
      if (party == ALICE) {
        uint64_t *min_index0 = new uint64_t[1];
        iopack->io->send_data(min_index, sizeof(uint64_t));
        iopack->io->recv_data(min_index0, sizeof(uint64_t));
        min_index[0] = (min_index[0] + min_index0[0]) & mask_x;
      } else {
        uint64_t *min_index0 = new uint64_t[1];
        iopack->io->send_data(min_index, sizeof(uint64_t));
        iopack->io->recv_data(min_index0, sizeof(uint64_t));
        min_index[0] = (min_index[0] + min_index0[0]) & mask_x;
      }
      clusters[min_index[0]].push_back(x[i]);
    }
    // update medians
    uint64_t *curr_cluster;
    for (int i = 0; i < K; i++) {
      curr_cluster = new uint64_t[clusters[i].size()];
      for (int j = 0; j < clusters[i].size(); j++) {
        curr_cluster[j] = clusters[i][j];
      }
      uint64_t *median_temp = new uint64_t[1];
      median(clusters[i].size(), curr_cluster, median_temp, bw_x);
      medians[i] = median_temp[0];
      delete[] curr_cluster;
      delete[] median_temp;
    }
  }

  for (int i = 0; i < K; i++) {
    y[i] = medians[i];
  }
}
