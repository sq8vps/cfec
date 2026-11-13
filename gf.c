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

#include "gf.h"
#include "gfprim.h"
#include <string.h>
#include <stdlib.h>

#define GF_EXTENSION_DEGREE_MAX 32


void GfPolyScale(const struct GaloisField *gf, 
    const DATA_TYPE *poly, size_t length, DATA_TYPE scalar, DATA_TYPE *out)
{
    for(size_t i = 0; i < length; i++)
    {
        out[i] = GfMul(gf, poly[i], scalar);
    }
}

void GfPolyScaleInPlace(const struct GaloisField *gf, 
    DATA_TYPE *poly, size_t length, DATA_TYPE scalar)
{
    GfPolyScale(gf, poly, length, scalar, poly);
}

size_t GfPolyAdd(const struct GaloisField *gf,
    const DATA_TYPE *poly1, size_t length1, 
    const DATA_TYPE *poly2, size_t length2, DATA_TYPE *out)
{
    size_t size = (length1 > length2) ? length1 : length2;
    size_t shorterSize = (length1 <= length2) ? length1 : length2;
    const DATA_TYPE *shorter = (length1 <= length2) ? poly1 : poly2;
    memcpy(out, (length1 > length2) ? poly1 : poly2, size * sizeof(*out));
    for(size_t i = 0; i < shorterSize; i++)
    {
        out[size - i - 1] = GfAdd(gf, out[size - i - 1], shorter[shorterSize - i - 1]);
    }
    return size;
}

size_t GfPolySub(const struct GaloisField *gf,
    const DATA_TYPE *poly1, size_t length1, 
    const DATA_TYPE *poly2, size_t length2, DATA_TYPE *out)
{
    if(length1 >= length2)
    {
        memcpy(out, poly1, length1 * sizeof(*out));
        for(size_t i = 0; i < length2; i++)
        {
            out[length1 - i - 1] = GfSub(gf, out[length1 - i - 1], poly2[length2 - i - 1]);
        }
    }
    else //length2 > length1
    {
        memset(out, 0, (length2 - length1) * sizeof(*out));
        memcpy(&out[length2 - length1], poly1, length1 * sizeof(*out));
        for(size_t i = 0; i < length2; i++)
        {
            out[i] = GfSub(gf, out[i], poly2[i]);
        }
    }
    return (length1 >= length2) ? length1 : length2;
}

void GfPolyMul(const struct GaloisField *gf,
    const DATA_TYPE *poly1, size_t length1, 
    const DATA_TYPE *poly2, size_t length2, DATA_TYPE *out)
{
    memset(out, 0, (length1 + length2 - 1) * sizeof(*out));
    for(size_t i = 0; i < length1; i++)
    {
        for(size_t j = 0; j < length2; j++)
        {
            out[i + j] = GfAdd(gf, out[i + j], GfMul(gf, poly1[i], poly2[j]));
        }
    }
}

DATA_TYPE GfPolyEval(const struct GaloisField *gf, 
    const DATA_TYPE *poly, size_t length, DATA_TYPE x)
{
    //this uses Horner's scheme to perform fast evaluation
    DATA_TYPE result = 0;
    result = poly[0];
    for(size_t i = 1; i < length; i++)
    {
        result = GfAdd(gf, GfMul(gf, result, x), poly[i]); 
    }
    return result;
}

DATA_TYPE *GfPolyDiv(const struct GaloisField *gf,
    const DATA_TYPE *poly1, size_t length1, 
    const DATA_TYPE *poly2, size_t length2, DATA_TYPE *out)
{
    DATA_TYPE normalizer = poly2[0];
    memcpy(out, poly1, length1 * sizeof(*out));
    for(size_t i = 0; i < (length1 - length2 + 1); i++)
    {
        out[i] = GfDiv(gf, out[i], normalizer);
        DATA_TYPE coeff = out[i];
        if(0 == coeff) 
            continue;

        for(size_t j = 1; j < length2; j++)
        {
            if(0 == poly2[j]) 
                continue;
            out[i + j] = GfSub(gf, out[i + j], GfMul(gf, poly2[j], coeff));
        }
    }
    return &(out[length1 - length2 + 1]);
}

void GfPolyInv(const struct GaloisField *gf, DATA_TYPE *poly, size_t length)
{
    for(size_t i = 0; i < (length / 2); i++)
    {
        DATA_TYPE tmp = poly[i];
        poly[i] = poly[length - i - 1];
        poly[length - i - 1] = tmp;
    }
}

void GfPolyNeg(const struct GaloisField *gf, DATA_TYPE *poly, size_t length)
{
    for(size_t i = 0; i < length; i++)
    {
        poly[i] = GfNeg(gf, poly[i]);
    }
}

/**
 * @brief Encode extension field coefficients into a single number
 * @param *coeffs Field element polynomial coefficients, highest degree first
 * @param p Field characteristics
 * @param n Field degree (number of coefficients)
 * @return Encoded number
 */
static DATA_TYPE GfEncodeExtension(const DATA_TYPE *coeffs, uint32_t p, uint32_t n)
{
    DATA_TYPE x = 0;
    DATA_TYPE mul = 1;
    for(size_t i = 0; i < n; i++)
    {
        x += coeffs[i] * mul;
        mul *= p;
    }
    return x;
}

/**
 * @brief Decode extension field coefficients from a single number
 * @param x Number to decode
 * @param *coeffs Field element polynomial coefficients, highest degree first
 * @param p Field characteristics
 * @param n Field degree (number of coefficients)
 */
static void GfDecodeExtension(DATA_TYPE x, DATA_TYPE *coeffs, uint32_t p, uint32_t n)
{
    for(size_t i = 0; i < n; i++)
    {
        coeffs[i] = x % p;
        x /= p;
    }
}

static DATA_TYPE GfAddInExtension(DATA_TYPE a, DATA_TYPE b, uint32_t p, uint32_t n)
{
    DATA_TYPE result = 0;
    DATA_TYPE mul = 1;

    for(size_t i = 0; i < n; i++)
    {
        DATA_TYPE ai = a % p;
        DATA_TYPE bi = b % p;

        result += ((uint64_t)(ai + bi) % p) * mul;

        mul *= p;
        a /= p;
        b /= p;
    }

    return result;
}

static DATA_TYPE GfNegInExtension(DATA_TYPE x, uint32_t p, uint32_t n)
{
    DATA_TYPE result = 0;
    DATA_TYPE mul = 1;

    for(size_t i = 0; i < n; i++)
    {
        DATA_TYPE xi = x % p;

        if(0 != xi)
            result += (p - xi) * mul;

        mul *= p;
        x /= p;
    }

    return result;    
}

static DATA_TYPE GfMulInExtension(DATA_TYPE a, DATA_TYPE b, const DATA_TYPE *modc, uint32_t p, uint32_t n)
{
    DATA_TYPE result[2 * GF_EXTENSION_DEGREE_MAX];
    DATA_TYPE ac[GF_EXTENSION_DEGREE_MAX];
    DATA_TYPE bc[GF_EXTENSION_DEGREE_MAX];

    memset(result, 0, sizeof(result));

    GfDecodeExtension(a, ac, p, n);
    GfDecodeExtension(b, bc, p, n);

    //calculate normal polynomial product first
    for(size_t i = 0; i < n; i++)
    {
        for(size_t j = 0; j < n; j++)
        {
            result[i + j] = ((uint64_t)result[i + j] + (uint64_t)ac[i] * (uint64_t)bc[j]) % p;
        }
    }

    //reduce the resulting product with primitive polynomial
    //this is the same as for GF(2^n), except here it's more complicated
    for(size_t i = (2 * n - 2); i >= n; i--)
    {
        DATA_TYPE c = result[i];

        if(0 != c)
        {
            for(size_t j = 0; j < n; j++)
            {
                DATA_TYPE d = ((uint64_t)c * (uint64_t)modc[j]) % p;
                result[i + j - n] = ((uint64_t)result[i + j - n] + (uint64_t)p - (uint64_t)d) % p;
            }
        }

        result[i] = 0;
    }

    return GfEncodeExtension(result, p, n);
}

static size_t GfFindPrimeFactors(uint64_t n, uint64_t factor[9])
{
    size_t count = 0;

    if(0 == (n % 2))
    {
        factor[count++] = 2;
        while(0 == (n % 2))
            n /= 2;
    }

    for(uint64_t i = 3; (i * i) <= n; i += 2)
    {
        if(0 == (n % i))
        {
            factor[count++] = i;
            while(0 == (n % i))
                n /= i;
        }
    }

    if(n > 1)
        factor[count++] = n;

    return count;
}

static uint64_t GfRand(uint64_t max)
{
    uint64_t random = (uint64_t)rand();
    random |= ((uint64_t)rand() << 15);
    random |= ((uint64_t)rand() << 30);
    random |= ((uint64_t)rand() << 45);
    random |= ((uint64_t)rand() << 60);

    return 2 + (random % (max - 1));
}

static uint64_t GfPowMod(uint64_t base, uint64_t exponent, uint64_t mod)
{
    uint64_t result = 1;
    while(exponent--)
    {
        result = ((uint64_t)result * (uint64_t)base) % (uint64_t)mod;
    }
    return result;
}

static DATA_TYPE GfPowModInExtension(DATA_TYPE base, DATA_TYPE exponent, const DATA_TYPE *modc, uint32_t p, uint32_t n)
{
    DATA_TYPE result = 1;

    while(exponent)
    {
        if(exponent & 1)
            result = GfMulInExtension(result, base, modc, p, n);

        exponent >>= 1;

        if(exponent)
            base = GfMulInExtension(base, base, modc, p, n);
    }

    return result;
}

/**
 * @brief Find or check primitive element
 * @param check If not zero, it is checked if it is a primitive element
 * @param order Field order (p^n)
 * @param *modc Generator polynomial for extension fields
 * @param p Field characteristics
 * @param n Field degree
 * @return If \a check is zero, returns a primitive element. 
 * If \a check is not zero and \a check is a primitive element, \a check is returned. Otherwise zero is returned.
 */
static uint64_t GfFindOrCheckPrimitiveElement(uint64_t check, uint64_t order, const DATA_TYPE *modc, uint32_t p, uint32_t n)
{
    //find prime factors of p-1
    uint64_t exponent[9];
    size_t exponentCount = GfFindPrimeFactors(order - 1, exponent);
    //calculate test exponents as floor((p-1)/prime factor)
    for(size_t i = 0; i < exponentCount; i++)
        exponent[i] = (order - 1) / exponent[i];
    
    while(1)
    {
        bool fail = false;
        uint64_t primitive = check ? check : GfRand(order - 1);
        for(size_t i = 0; i < exponentCount; i++)
        {
            if(1 != n)
            {
                if(1 == GfPowModInExtension(primitive, exponent[i], modc, p, n))
                {
                    fail = true;
                    break;
                }
            }
            else
            {
                if(1 == GfPowMod(primitive, exponent[i], order))
                {
                    fail = true;
                    break;
                }
            }
        }

        if(check)
            return fail ? 0 : check;
        if(!fail)
            return primitive;
    }
}

static bool GfIsPrime(uint64_t n)
{
    if(n <= 1)
        return false;

    if((2 == n) || (3 == n))
        return true;

    if((0 == (n % 2)) || (0 == (n % 3)))
        return false;
    
    for(uint64_t i = 5; (i * i) <= n; i += 6)
    {
        if((0 == (n % i)) || (0 == (n % (i + 2))))
            return false;
    }

    return true;
}

CFecStatus GfInit(uint32_t characteristic, uint32_t degree, uint32_t primitive, struct GaloisField **gfUser)
{
    struct GaloisField *gf = NULL;
    uint64_t order = 1;
    bool extension = !!((2 != characteristic) && (1 != degree)); //extension field other than 2^n?

    if(!GfIsPrime(characteristic) || (characteristic < 2) || (0 == degree) || (degree > GF_EXTENSION_DEGREE_MAX))
        return CFEC_BAD_PARAMETER;

    if(NULL == gfUser)
        return CFEC_BAD_PARAMETER;

    if(extension && (characteristic > GF_MAX_EXTENDED_CHARACTERISTICS))
        return CFEC_BAD_PARAMETER;

    for(uint32_t i = 0; i < degree; i++)
    {
        order *= characteristic;
        if(order > ((uint64_t)1 << (8 * sizeof(DATA_TYPE))))
            return CFEC_BAD_PARAMETER;
    }

    if(extension && ((order * order) <= order))
        return CFEC_NO_RESOURCES;

    if(extension)
    {
        gf = malloc(sizeof(*gf) + (3 * sizeof(*gf->lut.exp) * order) + (sizeof(*gf->lut.add) * order * order));
        if(NULL == gf)
            return CFEC_NO_RESOURCES;
        gf->lut.exp = (void*)((uintptr_t)gf + sizeof(*gf));
        gf->lut.log = (void*)((uintptr_t)gf->lut.exp + sizeof(*gf->lut.exp) * order);
        gf->lut.neg = (void*)((uintptr_t)gf->lut.log + sizeof(*gf->lut.log) * order);
        gf->lut.add = (void*)((uintptr_t)gf->lut.neg + sizeof(*gf->lut.neg) * order);
    }
    else
    {
        gf = malloc(sizeof(*gf) + (2 * sizeof(*gf->lut.exp) * order));
        if(NULL == gf)
            return CFEC_NO_RESOURCES;
        gf->lut.exp = (void*)((uintptr_t)gf + sizeof(*gf));
        gf->lut.log = (void*)((uintptr_t)gf->lut.exp + sizeof(*gf->lut.exp) * order);
    }

    gf->order = order;
    
    if(2 == characteristic) //GF(2^n)
    {
        uint64_t poly = Gf2PrimitivePoly[degree - 1];
        uint64_t x = 1;
        
        for(size_t i = 0; i < (order - 1); i++)
        {
            gf->lut.exp[i] = x;
            gf->lut.log[x] = i;

            x <<= 1;

            if(x & order)
                x ^= poly;
            x &= (order - 1);
        }

        gf->binary = true;
        gf->extension = false;
        gf->alpha = 2;
    }
    else //GF(p^n) where p>2
    {
        uint64_t generator = GfFindOrCheckPrimitiveElement(primitive, order, 
            extension ? GfPrimitivePoly[GfPrimitivePolyMap[characteristic]][degree - 2] : NULL, 
            characteristic, degree);
        if(0 == generator)
        {
            free(gf);
            return CFEC_BAD_PARAMETER;
        }
        gf->binary = false;
        gf->alpha = generator;
        gf->extension = extension;

        uint64_t x = 1;
        for(size_t i = 0; i < (order - 1); i++)
        {
            gf->lut.exp[i] = x;
            gf->lut.log[x] = i;

            if(extension)
                x = GfMulInExtension(generator, x, 
                    GfPrimitivePoly[GfPrimitivePolyMap[characteristic]][degree - 2], characteristic, degree);
            else
                x = (x * generator) % (uint64_t)order;
        }

        if(extension)
        {
            for(uint64_t i = 0; i < order; i++)
            {
                for(uint64_t j = 0; j < order; j++)
                {
                    gf->lut.add[i * order + j] = GfAddInExtension(i, j, characteristic, degree);
                }

                gf->lut.neg[i] = GfNegInExtension(i, characteristic, degree);
            }
        }
    }

    *gfUser = gf;
    return CFEC_OK;
}

void GfDeinit(struct GaloisField *gf)
{
    free(gf);
}