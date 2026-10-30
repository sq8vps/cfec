#ifndef GF_H_
#define GF_H_

#include <stdint.h>
#include <stddef.h>

#define likely(x) (x)
#define unlikely(x) (x)

/**
 * @brief Galois Field library structure 
 */
struct GaloisField
{
    size_t order; /**< Field order */
    size_t wrap; /**< Value used as for masking and wrapping */
    bool binary; /**< Is order a power of 2? */
    uint32_t *expTable; /**< Precalculated table of exp(x) */
    uint32_t *logTable; /**< Precalculated table of log(x) */
};

/**
 * @brief Add in Galois field
 * @param *gf Galois field instance
 * @param x Term 1
 * @param y Term 2
 * @return Sum
 */
static inline uint32_t GfAdd(const struct GaloisField *gf, uint32_t x, uint32_t y)
{
    if(gf->binary)
        return (x ^ y) & gf->wrap;
    else
        return (x + y) % gf->wrap;
}

/**
 * @brief Subtract in Galois field
 * @param *gf Galois field instance
 * @param x Minuend
 * @param y Subtrahend
 * @return Difference
 */
static inline uint32_t GfSub(const struct GaloisField *gf, uint32_t x, uint32_t y)
{
    if(gf->binary)
        return (x ^ y) & gf->wrap;
    else
        return (x - y) % gf->wrap;
}

/**
 * @brief Multiply in Galois field
 * @param *gf Galois field instance
 * @param x Multiplicand
 * @param y Multiplier
 * @return Multiplication result
 * @attention 0 is returned when \a x or \a y is bigger or equal to GF order
 */
static inline uint32_t GfMul(const struct GaloisField *gf, uint32_t x, uint32_t y)
{
    if((0 == x) || (0 == y))
        return 0;

    if(unlikely((x >= gf->order) || (y >= gf->order)))
        return 0;

    //use b^(log(x)+log(y))=b^log(x*y)=x*y identity
    if(gf->binary)
        return gf->expTable[(gf->logTable[x] + gf->logTable[y]) & gf->wrap];
    else
        return gf->expTable[(gf->logTable[x] + gf->logTable[y]) % gf->wrap];
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
static inline uint32_t GfDiv(const struct GaloisField *gf, uint32_t dividend, uint32_t divisor)
{
    if(unlikely(0 == divisor)) 
        return 0;

    if(unlikely((dividend >= gf->order) || (divisor >= gf->order)))
        return 0;
        
    if(0 == dividend)
        return 0;

    //use b^(log(x)-log(y))=b^(log(x/y))=x/y identity
    if(gf->binary)
        return gf->expTable[(gf->order - 1 + gf->logTable[dividend] - gf->logTable[divisor]) & gf->wrap];
    else
        return gf->expTable[(gf->order + gf->logTable[dividend] - gf->logTable[divisor]) % gf->wrap];
}

/**
 * @brief Exponentiate in Galois field
 * @param *gf Galois field instance
 * @param x Base
 * @param exponent Exponent
 * @return Result
 */
static inline uint32_t GfPow(const struct GaloisField *gf, uint32_t x, uint32_t exponent)
{
    if(unlikely((x >= gf->order) || (exponent >= gf->order)))
        return 0;
        
    if(0 == exponent)
        return 1;

    //use b^(a*log(x))=b^(log(x^a))=x^a identity
    if(gf->binary)
        return gf->expTable[(exponent * GfLog[x]) & gf->wrap];
    else
        return gf->expTable[(exponent * GfLog[x]) % gf->wrap];
}

/**
 * @brief Invert in Galois field
 * @param *gf Galois field instance
 * @param x Number to calculate the inverse of
 * @return 1/x
 */
static inline uint32_t GfInv(const struct GaloisField *gf, uint32_t x)
{
    if(unlikely(x >= gf->order))
        return 0;
    
    return gf->expTable[gf->wrap - gf->logTable[x]];
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
    const uint32_t *poly, size_t length, uint32_t scalar, uint32_t *out);

/**
 * @brief Multiply polynomial by a scalar in-place
 * @param *gf Galois field instance
 * @param *poly Polynomial
 * @param length Polynomial length (degree + 1)
 * @param scalar Scalar multiplicand
 */
void GfPolyScaleInPlace(const struct GaloisField *gf, 
    uint32_t *poly, size_t length, uint32_t scalar);

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
uint32_t GfPolyAdd(const struct GaloisField *gf,
    const uint32_t *poly1, size_t length1, 
    const uint32_t *poly2, size_t length2, uint32_t *out);

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
    const uint32_t *p1, size_t length1, 
    const uint32_t *poly2, size_t length2, uint32_t *out);

/**
 * @brief Evaluate the polynomial at given x
 * @param *gf Galois field instance
 * @param *poly Polynomial to evaulate
 * @param length Polynomial length (degree + 1)
 * @param x Value to evaluate the polynomial at
 * @return Evaluated value
 */
uint32_t GfPolyEval(const struct GaloisField *gf, 
    const uint32_t *poly, size_t length, uint32_t x);

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

uint32_t *GfPolyDiv(const struct GaloisField *gf,
    const uint32_t *poly1, size_t length1, 
    const uint32_t *poly2, size_t length2, uint32_t *out);

/**
 * @brief Reverse the order of elements in a polynomial (in-place)
 * @param *gf Galois field instance
 * @param *poly Polynomial input/output
 * @param length Polynomial length (degree + 1)
 */
void GfPolyInv(const struct GaloisField *gf, const uint32_t *poly, size_t length);

#endif