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

/*
Copyright 2023-2026 Piotr Wilkoń
This file is part of LwFEC.

LwFEC is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

LwFEC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LwFEC.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include "cfec.h"

#ifndef DATA_TYPE
#define DATA_TYPE uint8_t
#error RS: Data type not defined
#endif

struct GaloisField;

/**
 * @brief Reed-Solomon module configuration and operation structure
*/
struct CFecRS
{
    struct GaloisField *gf; /**< Galois field instance */
    size_t N; /**< Block size */
    size_t T; /**< Number of redundancy/parity bytes */
    size_t fcr; /**< First consecutive root index */
    DATA_TYPE generator[]; /**< Generator polynomial */
};

/**
 * @brief Decode message using Reed-Solomon FEC
 * @param *rs RS coder/decoder instance
 * @param *data Input/output buffer
 * @param size Data size
 * @param *fixed Output number of bytes corrected. Can be set to \a NULL if not needed.
 * @return Status code
 * 
 * This function takes input buffer with \a size data bytes and \a T parity bytes
 * Then it moves parity bytes to the end and fills everything inbetween with zeros.
 * Next the in-place decoding is performed.
 */
CFecStatus RsDecode(struct CFecRS *rs, DATA_TYPE *data, size_t size, size_t *fixed);

/**
 * @brief Encode message using Reed-Solomon FEC
 * @param *rs RS coder/decoder instance
 * @param *data Input/output buffer. Must have at least \a T more free elements than \a size.
 * @param size Data size
 * @return Status code
 */
CFecStatus RsEncode(struct CFecRS *rs, DATA_TYPE *data, size_t size);

/**
 * @brief Initialize Reed-Solomon coder/decoder
 * @param N Block size = dictionary size
 * @param T Number of parity check elements
 * @param fcr First consecutive root index
 * @param **rs Pointer to where to return the allocated Reed-Solomon instance
 * @return Status code
*/
CFecStatus RsInit(size_t N, size_t T, size_t fcr, struct CFecRS **rs);
