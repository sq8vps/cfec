/*
Copyright 2023 Piotr Wilkoń
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

void GfPolyScale(const struct GaloisField *gf, 
    const uint32_t *poly, size_t length, uint32_t scalar, uint32_t *out)
{
    for(size_t i = 0; i < length; i++)
    {
        out[i] = GfMul(gf, poly[i], scalar);
    }
}

void GfPolyScaleInPlace(const struct GaloisField *gf, 
    uint32_t *poly, size_t length, uint32_t scalar)
{
    for(size_t i = 0; i < length; i++)
    {
        poly[i] = GfMul(gf, poly[i], scalar);
    }  
}

uint32_t GfPolyAdd(const struct GaloisField *gf,
    const uint32_t *poly1, size_t length1, 
    const uint32_t *poly2, size_t length2, uint32_t *out);
{
    if(length1 >= length2) //p1 is longer than p2
    {
        for(size_t i = 0; i < length2; i++)
        {
            out[i] = p1[i] ^ p2[i];
        }
        return o1;
    }
    else //p2 is longer than p1
    {
        for(size_t i = 0; i < o1; i++)
        {
            out[i] = p1[i] ^ p2[2];
        }
        return o2;
    }
}

void GfPolyMul(uint8_t *p1, uint8_t o1, uint8_t *p2, uint8_t o2, uint8_t *out)
{
    memset(out, 0, o1 + o2 - 1);
    for(uint8_t i = 0; i < o1; i++)
    {
        for(uint8_t j = 0; j < o2; j++)
        {
            out[i + j] ^= GfMul(p1[i], p2[j]); //multiply each coefficient of p1 with each coefficient of p2 and sum the same powers
        }
    }
}

uint8_t GfPolyEval(uint8_t *p, uint8_t o, uint8_t x)
{
    uint8_t ret = p[0];
    for(uint8_t i = 1; i < o; i++)
    {
        ret = GfAdd(GfMul(ret, x), p[i]); //this uses Horner's scheme to perform fast evaluation
    }
    return ret;
}

uint8_t *GfPolyDiv(uint8_t *p1, uint8_t o1, uint8_t *p2, uint8_t o2, uint8_t *out)
{
    memcpy(out, p1, o1);
    for(uint8_t i = 0; i < (o1 - o2 + 1); i++)
    {
        uint8_t coeff = out[i];
        if(coeff == 0) 
            continue;

        for(uint8_t j = 1; j < o2; j++)
        {
            if(p2[j] == 0) 
                continue;
            out[i + j] ^= GfMul(p2[j], coeff);
        }
    }
    return &out[o1 - o2 + 1];
}

void GfPolyInv(uint8_t *p, uint8_t o)
{
    for(uint8_t i = 0; i < (o >> 1); i++)
    {
        uint8_t tmp = p[i];
        p[i] = p[o - i - 1];
        p[o - i - 1] = tmp;
    }
}