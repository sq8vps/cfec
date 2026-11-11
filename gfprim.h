/*
Copyright 2026 Piotr Wilkoń
This file is part of CFEC.

CFEC is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

CFEC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CFEC.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GFPRIM_H_
#define GFPRIM_H_

#include <stdint.h>
#include <stddef.h>

#ifndef DATA_TYPE
#define DATA_TYPE uint8_t
#error GF: Data type must be defined
#endif

/**
 * @brief Primitive polynomials for GF(2^k)
 */
static const uint64_t Gf2PrimitivePoly[32] = {
    0x3,
    0x7,
    0xB,
    0x13,
    0x25,
    0x43,
    0x83,
    0x11D,
    0x211,
    0x409,
    0x805,
    0x1053,
    0x201B,
    0x402B,
    0x8003,
    0x1002D,
    0x20009,
    0x40081,
    0x80027,
    0x100009,
    0x200005,
    0x400003,
    0x800021,
    0x100001B,
    0x2000009,
    0x4000047,
    0x8000027,
    0x10000009,
    0x20000005,
    0x40000053,
    0x80000009,
    0x190000003,
};

/**
 * @brief Table of primitive polynomials for different field extensions
 * @note The 3rd argument is the (degree - 2), because no polynomials are stored for degrees 0 and 1
 */
static const DATA_TYPE GfPrimitivePoly[][19][21] = {
    //GF(3^n)
    {
    {2, 1, 1},
    {1, 0, 2, 1},
    {2, 0, 0, 1, 1},
    {1, 0, 0, 0, 2, 1},
    {2, 0, 0, 0, 0, 1, 1},
    {1, 0, 0, 0, 0, 1, 2, 1},
    {2, 0, 0, 0, 0, 1, 0, 0, 1},
    {1, 0, 0, 0, 0, 0, 2, 1, 0, 1},
    {2, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1},
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1},
    {2, 0, 0, 1, 2, 2, 0, 0, 2, 1, 2, 1, 1},
    {1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 2, 1},
    {2, 0, 0, 2, 2, 0, 0, 2, 1, 1, 1, 1, 1, 0, 1},
    {1, 0, 1, 1, 2, 1, 2, 2, 2, 1, 0, 0, 2, 1, 2, 1},
    {2, 1, 2, 2, 0, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2, 2, 1},
    {1, 2, 2, 0, 0, 2, 0, 2, 1, 2, 2, 2, 1, 0, 2, 2, 0, 1},
    {2, 1, 0, 0, 0, 2, 1, 0, 0, 1, 0, 2, 0, 0, 0, 1, 1, 2, 1},
    {1, 1, 2, 0, 2, 2, 1, 0, 2, 0, 2, 0, 1, 2, 1, 1, 0, 2, 2, 1},
    {2, 0, 0, 2, 0, 1, 0, 2, 2, 1, 1, 0, 0, 0, 0, 2, 0, 2, 1, 2, 1}},
    //GF(5^n)
    {
    {3, 2, 1},
    {2, 0, 3, 1},
    {3, 3, 0, 2, 1},
    {3, 2, 3, 0, 3, 1},
    {2, 4, 1, 3, 4, 4, 1},
    {3, 3, 1, 2, 0, 2, 2, 1},
    {2, 4, 3, 1, 1, 1, 4, 4, 1},
    {2, 4, 0, 1, 2, 2, 2, 2, 2, 1},
    {3, 1, 3, 1, 1, 0, 1, 0, 2, 0, 1},
    {3, 1, 0, 2, 1, 2, 4, 1, 2, 3, 4, 1},
    {2, 3, 0, 0, 2, 2, 4, 0, 0, 1, 2, 1, 1},
    {3, 2, 1, 1, 0, 1, 4, 3, 3, 2, 3, 0, 2, 1}},
    //GF(7^n)
    {
    {5, 3, 1},
    {4, 5, 3, 1},
    {3, 5, 6, 0, 1},
    {4, 5, 3, 0, 5, 1},
    {5, 2, 0, 4, 2, 3, 1},
    {4, 3, 1, 2, 1, 1, 0, 1},
    {5, 5, 2, 4, 2, 0, 4, 0, 1},
    {4, 0, 1, 5, 1, 3, 1, 2, 0, 1},
    {3, 6, 0, 4, 6, 0, 5, 0, 4, 0, 1},
    {2, 3, 3, 3, 5, 4, 4, 0, 1, 3, 1, 1}},
    //GF(11^n)
    {
    {7, 1, 1},
    {3, 0, 1, 1},
    {8, 0, 0, 1, 1},
    {3, 0, 0, 1, 1, 1},
    {8, 0, 0, 0, 1, 5, 1},
    {4, 0, 0, 0, 0, 0, 1, 1},
    {7, 0, 0, 0, 0, 0, 1, 1, 1},
    {3, 0, 0, 0, 0, 0, 0, 0, 2, 1}},
    //GF(13^n)
    {
    {2, 1, 1},
    {2, 0, 1, 1},
    {6, 0, 1, 1, 1},
    {7, 0, 0, 0, 2, 1},
    {11, 0, 0, 0, 1, 5, 1},
    {6, 0, 0, 0, 0, 0, 5, 1},
    {11, 0, 0, 0, 0, 0, 1, 1, 1}},
    //GF(17^n)
    {
    {3, 1, 1},
    {7, 0, 1, 1},
    {5, 0, 0, 1, 1},
    {5, 0, 0, 0, 1, 1},
    {3, 0, 0, 0, 0, 1, 1},
    {7, 0, 0, 0, 0, 0, 1, 1}}
};

/**
 * @brief Mapping of the extension field characteristics to the appropriate primitive poly set
 */
static const size_t GfPrimitivePolyMap[] = {0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 3, 0, 4, 0, 0, 0, 5};

/**
 * @brief Highest supported extended field characteristics
 */
#define GF_MAX_EXTENDED_CHARACTERISTICS 17

#endif