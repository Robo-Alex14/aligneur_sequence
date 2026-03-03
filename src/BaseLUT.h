#pragma once
#include <cstdint>

// ─────────────────────────────────────────────────────────────────────────────
// BASE_LUT: 256-entry lookup table mapping every possible byte value to a
// 2-bit base code or INVALID.
//
//   A / a  →  0
//   C / c  →  1
//   G / g  →  2
//   T / t  →  3
//   anything else (N, n, \r, space, …)  →  INVALID (0xFF)
//
// Using a direct array lookup (O(1), branch-free) instead of a switch/toupper
// chain means that ANY unexpected character — Windows \r, lowercase, spaces,
// IUPAC ambiguity codes, etc. — is automatically treated as invalid without
// any special-casing.
//
// All modules that need to encode or compare bases #include this header and
// call BASE_LUT[c] where c is an unsigned char.
// ─────────────────────────────────────────────────────────────────────────────

static const uint8_t INVALID = 0xFF;

// Defined in BaseLUT.cpp; declared extern so every translation unit that
// includes this header shares the same table.
extern const uint8_t BASE_LUT[256];
