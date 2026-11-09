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

#ifndef GF_H_
#define GF_H_

#include <stdint.h>
#include <stddef.h>
#include "cfec.h"

#ifndef likely
#define likely(x) (x)
#endif

#ifndef unlikely
#define unlikely(x) (x)
#endif

#ifndef DATA_TYPE
#define DATA_TYPE uint8_t
#error GF: Data type not defined
#endif

/**
 * @brief Galois Field library structure 
 */
struct GaloisField
{
    uint64_t order; /**< Field order */
    DATA_TYPE alpha; /**< Field generator */
    bool binary; /**< Is order a power of 2? */
    DATA_TYPE *expTable; /**< Precalculated table of exp(x) */
    DATA_TYPE *logTable; /**< Precalculated table of log(x) */
    char data[]; /**< Precalculated tabled will be stored starting from here */
};

/**
 * @brief Add in Galois field
 * @param *gf Galois field instance
 * @param x Term 1
 * @param y Term 2
 * @return Sum
 */
static inline DATA_TYPE GfAdd(const struct GaloisField *gf, DATA_TYPE x, DATA_TYPE y)
{
    if(gf->binary)
        return (x ^ y) & (gf->order - 1);
    else
        return (x + y) % gf->order;
}

/**
 * @brief Subtract in Galois field
 * @param *gf Galois field instance
 * @param x Minuend
 * @param y Subtrahend
 * @return Difference
 */
static inline DATA_TYPE GfSub(const struct GaloisField *gf, DATA_TYPE x, DATA_TYPE y)
{
    if(gf->binary)
        return (x ^ y) & (gf->order - 1);
    else
        return (gf->order + x - y) % gf->order;
}

/**
 * @brief Negate number in Galois field
 * @param *gf Galois field instance
 * @param x Number
 * @return Number negated
 */
static inline DATA_TYPE GfNeg(const struct GaloisField *gf, DATA_TYPE x)
{
    if(gf->binary)
        return x;
    else
    {
        if(x >= gf->order)
            x %= gf->order;
        return (gf->order - x) % gf->order;
    }
}

/**
 * @brief Multiply in Galois field
 * @param *gf Galois field instance
 * @param x Multiplicand
 * @param y Multiplier
 * @return Multiplication result
 * @attention 0 is returned when \a x or \a y is bigger or equal to GF order
 */
static inline DATA_TYPE GfMul(const struct GaloisField *gf, DATA_TYPE x, DATA_TYPE y)
{
    if((0 == x) || (0 == y))
        return 0;

    if(x >= gf->order)
        x %= gf->order;

    if(y >= gf->order)
        y %= gf->order;

    //use b^(log(x)+log(y))=b^log(x*y)=x*y identity
    return gf->expTable[(gf->logTable[x] + gf->logTable[y]) % (gf->order - 1)];
}

/**
 * @brief Divide in Galois field
 * @param *gf Galois field instance
 * @param dividend Dividend
 * @param divisor Divisor
 * @return Division result
 * @attention 0 is returned when \a x or \a y is bigger or equal to GF order
 * @attention 0 is returned when dividing by zero
 */
static inline DATA_TYPE GfDiv(const struct GaloisField *gf, DATA_TYPE dividend, DATA_TYPE divisor)
{
    if(unlikely(0 == divisor)) 
        return 0;

    if(dividend >= gf->order)
        dividend %= gf->order;
    
    if(divisor >= gf->order)
        divisor %= gf->order;
        
    if(0 == dividend)
        return 0;

    //use b^(log(x)-log(y))=b^(log(x/y))=x/y identity
    return gf->expTable[((gf->order - 1) + gf->logTable[dividend] - gf->logTable[divisor]) % (gf->order - 1)];
}

/**
 * @brief Exponentiate in Galois field
 * @param *gf Galois field instance
 * @param x Base
 * @param exponent Exponent
 * @return Result
 */
static inline DATA_TYPE GfPow(const struct GaloisField *gf, DATA_TYPE x, DATA_TYPE exponent)
{
    if(x >= gf->order) 
        x %= gf->order;
    
    if(exponent >= gf->order)
        exponent %= gf->order;
        
    if(0 == exponent)
        return 1;

    if(0 == x)
        return 0;

    //use b^(a*log(x))=b^(log(x^a))=x^a identity
    return gf->expTable[(exponent * gf->logTable[x]) % (gf->order - 1)];
}

/**
 * @brief Exponentiate field generator in Galois field
 * @param *gf Galois field instance
 * @param exponent Exponent
 * @return Generator number (alpha) to the power of \a exponent
 */
static inline DATA_TYPE GfPowAlpha(const struct GaloisField *gf, DATA_TYPE exponent)
{
    if(exponent >= gf->order)
        exponent %= gf->order;
        
    if(0 == exponent)
        return 1;

    return gf->expTable[(exponent * gf->logTable[gf->alpha]) % (gf->order - 1)];
}

/**
 * @brief Invert in Galois field
 * @param *gf Galois field instance
 * @param x Number to calculate the inverse of
 * @return 1/x
 */
static inline DATA_TYPE GfInv(const struct GaloisField *gf, DATA_TYPE x)
{
    if(x >= gf->order)
        x %= gf->order;
    
    if(unlikely(0 == x))
        return 0;
    
    return gf->expTable[(gf->order - 1) - gf->logTable[x]];
}

/**
 * @brief Multiply polynomial by a scalar
 * @param *gf Galois field instance
 * @param *poly Input polynomial
 * @param length Polynomial length (degree + 1)
 * @param scalar Scalar multiplicand
 * @param out Ouput polynomial. Has the same degree as input polynomial
 */
void GfPolyScale(const struct GaloisField *gf, 
    const DATA_TYPE *poly, size_t length, DATA_TYPE scalar, DATA_TYPE *out);

/**
 * @brief Multiply polynomial by a scalar in-place
 * @param *gf Galois field instance
 * @param *poly Polynomial
 * @param length Polynomial length (degree + 1)
 * @param scalar Scalar multiplicand
 */
void GfPolyScaleInPlace(const struct GaloisField *gf, 
    DATA_TYPE *poly, size_t length, DATA_TYPE scalar);

/**
 * @brief Add two polynomials
 * @param *gf Galois field instance
 * @param *poly1 1st polynomial
 * @param length1 1st polynomial length (degree + 1)
 * @param *poly2 2nd polynomial
 * @param length2 2nd polynomial length (degree + 1)
 * @param *out Output polynomial buffer
 * @return Output polynomial length
 */
size_t GfPolyAdd(const struct GaloisField *gf,
    const DATA_TYPE *poly1, size_t length1, 
    const DATA_TYPE *poly2, size_t length2, DATA_TYPE *out);

/**
 * @brief Subtract two polynomials
 * @param *gf Galois field instance
 * @param *poly1 Minuend polynomial
 * @param length1 Minuend polynomial length (degree + 1)
 * @param *poly2 Subtrahend polynomial
 * @param length2 Subtrahend polynomial length (degree + 1)
 * @param *out Output polynomial buffer
 * @return Output polynomial length
 */
size_t GfPolySub(const struct GaloisField *gf,
    const DATA_TYPE *poly1, size_t length1, 
    const DATA_TYPE *poly2, size_t length2, DATA_TYPE *out);

/**
 * @brief Multiply two polynomials
 * @param *gf Galois field instance
 * @param *poly1 1st polynomial
 * @param length1 1st polynomial length (degree + 1)
 * @param *poly2 2nd polynomial
 * @param length2 2nd polynomial length (degree + 1)
 * @param out Output polynomial. Has a length of \a length1 + \a length2 - 1
 * @warning Output buffer must be separate from input buffers
 */
void GfPolyMul(const struct GaloisField *gf,
    const DATA_TYPE *poly1, size_t length1, 
    const DATA_TYPE *poly2, size_t length2, DATA_TYPE *out);

/**
 * @brief Evaluate the polynomial at given x
 * @param *gf Galois field instance
 * @param *poly Polynomial to evaulate
 * @param length Polynomial length (degree + 1)
 * @param x Value to evaluate the polynomial at
 * @return Evaluated value
 */
DATA_TYPE GfPolyEval(const struct GaloisField *gf, 
    const DATA_TYPE *poly, size_t length, DATA_TYPE x);

/**
 * @brief Divide two polynomials in Galois field
 * @param *gf Galois field instance
 * @param *poly1 Divident polynomial
 * @param length1 Divident polynomial length (degree + 1)
 * @param *poly2 Divisor polynomial
 * @param length2 Divisor polynomial length (degree + 1)
 * @param *out Output polynomial quotient and remainder. Size \a length1 + \a length2 - 1. Must be preallocated.
 * @warning This function works on polynomials orderder highest-degree-term-first
 * @return Pointer to the first element of the remainder
 */

DATA_TYPE *GfPolyDiv(const struct GaloisField *gf,
    const DATA_TYPE *poly1, size_t length1, 
    const DATA_TYPE *poly2, size_t length2, DATA_TYPE *out);

/**
 * @brief Reverse the order of elements in a polynomial (in-place)
 * @param *gf Galois field instance
 * @param *poly Polynomial input/output
 * @param length Polynomial length (degree + 1)
 */
void GfPolyInv(const struct GaloisField *gf, DATA_TYPE *poly, size_t length);

/**
 * @brief Negate coefficients of the polynomial (in-place)
 * @param *gf Galois field instance
 * @param *poly Polynomial input/output
 * @param length Polynomial length (degree + 1)
 */
void GfPolyNeg(const struct GaloisField *gf, DATA_TYPE *poly, size_t length);

/**
 * @brief Initialize Galois Field library instance
 * @param order GF order, either a prime number \a p, or \a 2^k where 2<=k<=32
 * @param **gf Pointer to where to return the allocated Galois Field structure
 * @return Status code
 */
CFecStatus GfInit(uint64_t order, struct GaloisField **gf);

/**
 * @brief Deinitialize and remove Galois Field instance
 * @param *gf GF instance obtained from \ref GfInit()
 */
void GfDeinit(struct GaloisField *gf);

#endif