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
#include <string.h>
#include <stdlib.h>

static const uint64_t GfPrimitivePoly[32] = {
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

    return (random + 2) % max;
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

static __uint64_t GfFindPrimitiveElement(uint64_t order)
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
        uint64_t primitive = GfRand(order - 1);
        for(size_t i = 0; i < exponentCount; i++)
        {
            if(1 == GfPowMod(primitive, exponent[i], order))
            {
                fail = true;
                break;
            }
        }
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

static size_t GfOnePosition(uint64_t x)
{
    size_t p = 0;
    while(0 == (x & 1))
    {
        ++p;
        x >>= 1;
    }
    return p;
}

CFecStatus GfInit(uint64_t order, struct GaloisField **gfUser)
{
    struct GaloisField *gf = NULL;

    if(order <= 1)
        return CFEC_BAD_PARAMETER;

    if(NULL == gfUser)
        return CFEC_BAD_PARAMETER;

    gf = malloc(sizeof(*gf) + (2 * sizeof(*gf->expTable) * order));
    if(NULL == gf)
        return CFEC_NO_RESOURCES;

    gf->expTable = (void*)((uintptr_t)gf + sizeof(*gf));
    gf->logTable = (void*)((uintptr_t)gf->expTable + sizeof(*gf->expTable) * order);

    gf->order = order;
    
    if(0 == (order & (order - 1))) //is a power of 2?
    {
        if(order > ((uint64_t)1 << (8 * sizeof(DATA_TYPE))))
        {
            free(gf);
            return CFEC_BAD_PARAMETER;
        }

        uint64_t poly = GfPrimitivePoly[GfOnePosition(order) - 1];
        uint64_t x = 1;
        
        for(size_t i = 0; i < (order - 1); i++)
        {
            gf->expTable[i] = x;
            gf->logTable[x] = i;

            x <<= 1;

            if(x & order)
                x ^= poly;
            x &= (order - 1);
        }

        gf->binary = true;
        gf->alpha = 2;
    }
    else //not a power of 2
    {
        if(!GfIsPrime(order))
        {
            free(gf);
            return CFEC_BAD_PARAMETER;
        }

        uint64_t generator = GfFindPrimitiveElement(order);
        uint64_t x = 1;
        
        for(size_t i = 0; i < (order - 1); i++)
        {
            gf->expTable[i] = x;
            gf->logTable[x] = i;

            x = (x * generator) % (uint64_t)order;
        }

        gf->binary = false;
        gf->alpha = generator;
    }

    *gfUser = gf;
    return CFEC_OK;
}

void GfDeinit(struct GaloisField *gf)
{
    free(gf);
}