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

#include "rs.h"
#include "gf.h"
#include <string.h>
#include <stdlib.h>

//#define RS_USE_ALTERNATIVE_BM //use alternative Berlekamp-Massey implementation. Seems to be a bit faster

/**
 * @brief Calculates message syndromes
 * @param *rs RS instance
 * @param data Input block (length = N)
 * @param size Block size = N
 * @param out Output syndromes (length = T)
 * @note The output polynomial is the lowest degree
 */
static void RsSyndromes(const struct CFecRS *rs, const DATA_TYPE *data, size_t size, DATA_TYPE *out)
{
    //The syndrome polynomial is just the message evaluated at alpha^(i+fcr)
    //This produces a polynomial with the lowest degree first
    for(size_t i = 0; i < rs->T; i++)
    {
        out[i] = GfPolyEval(rs->gf, data, size, GfPowAlpha(rs->gf, i + rs->fcr));
    }
}

/**
 * @brief Obtain the list of error positions
 * @param *rs RS instance
 * @param *locator Error locator polynomial from \ref RsErrorLocator
 * @param locatorSize Error locator polynomial length
 * @param out List of erroneous positions (length = locatorSize - 1)
 * @return \a True on success, else the \a out buffer must be invalidated and the block is uncorrectable
 */
static bool RsErrorPositions(const struct CFecRS *rs, const DATA_TYPE *locator, size_t locatorSize, DATA_TYPE *out)
{
    //This is just finding the zeros of the error locator polynomial
    //These zeros are in form of alpha^i, where i denotes the error position
    //This can be done efficiently with Chien search, but it would require a buffer to store the coefficients
    //between iterations. For now let's stick to a brute-force approach with GfPolyEval() 
    size_t pos = 0;
    for(size_t i = 0; i < rs->N; i++)
    {
        if(GfPolyEval(rs->gf, locator, locatorSize, GfPowAlpha(rs->gf, i)) == 0) //if alpha^i is the root of the locator polynomial, it determines the error position
        {
            if(pos < (locatorSize - 1))
                out[pos] = rs->N - i - 1;
            else
                break;
            pos++;
        }
    }

    if(pos != (locatorSize - 1))
        return CFEC_UNCORRECTABLE;

    return CFEC_OK;
}


/**
 * @brief Calculate the error locator polynomial
 * @param *rs RS instance
 * @param *syndromes Syndrome polynomial obtained from \ref RsSyndromes()
 * @param *out Output error locator buffer
 * @param outSize Error locator polynomial buffer length
 * @return Status code
 */
static CFecStatus RsErrorLocator(struct CFecRS *rs, const DATA_TYPE *syndromes, DATA_TYPE *out, size_t *outSize)
{
    CFecStatus status = CFEC_OK;
    /*
     * The error locator polynomial is calculated using Berlekamp-Massey algorithm.
     * Two implementations are written here:
     * 1. directly adapted from Wikipedia (https://en.wikipedia.org/wiki/Berlekamp%E2%80%93Massey_algorithm)
     * 2. adapted/ported from Python from Wikiversity "Reed-Solomon codes for coders" 
     * (https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders#Error_correction)
     * Both implementations work just fine.
    */
#ifndef RS_USE_ALTERNATIVE_BM

    DATA_TYPE *B = calloc(rs->T + 1, sizeof(*B)); 
    DATA_TYPE *C = calloc(rs->T + 1, sizeof(*C)); //C(x) - candidate for error locator Lambda(x)
    DATA_TYPE *T = calloc(rs->T + 1, sizeof(*T));
    DATA_TYPE *T2 = calloc(rs->T + 1, sizeof(*T2));

    if(!B || !C || !T || !T2)
    {
        status = CFEC_NO_RESOURCES;
        goto leave;
    }

    size_t L = 0; //number of assumed errors
    size_t m = 1; //number of iterations since L, B and b were updated
    size_t b = 1; //last discrepancy delta

    //initialize B and C polynomials with the constant term
    B[0] = 1;
    C[0] = 1;

    for(size_t i = 0; i < rs->T; i++)
    {
        DATA_TYPE d = syndromes[i];
        for(size_t j = 1; j <= L; j++)
        {
            d = GfAdd(rs->gf, d, GfMul(rs->gf, C[j], syndromes[i - j])); //calculate discrepancy delta
        }
        if(d == 0) //no error
        {
            m++;
        }
        else if((L << 1) <= i)
        {
            memcpy(T, C, rs->T * sizeof(*T)); //store C polynomial in T
            //in general, C(x)=C(x)-(d/b)*B(x)*x^m
            //first B(x)=B(x)*x^m
            //here we store polynomials as the lowest degree term first
            //multipyling it by x^m shifts the polynomial coefficient by m places right
            //so swap places starting from the last element and fill first m places with zeros
            for(size_t j = 0; j < (rs->T - m); j++)
            {
                B[rs->T - j - 1] = B[rs->T - j - m - 1];
            }
            for(size_t j = 0; j < m; j++)
                B[j] = 0;

            //then B(x)*d/b
            GfPolyScale(rs->gf, B, rs->T, GfDiv(rs->gf, d, b), B);

            //then C(x)=C(x)-B(x)
            GfPolySub(rs->gf, T, rs->T, B, rs->T, C);

            //store T polynomial in B (previous C to B)
            memcpy(B, T, rs->T * sizeof(*B));

            L = i + 1 - L;
            b = d;
            m = 1;
        }
        else
        {
            //the same as above, but save B in T2 first
            memcpy(T, C, rs->T * sizeof(*T));
            memcpy(T2, B, rs->T * sizeof(*T));
            for(size_t j = 0; j < (rs->T - m); j++)
            {
                B[rs->T - j - 1] = B[rs->T - j - m - 1];
            }
            for(size_t j = 0; j < m; j++)
                B[j] = 0;
            GfPolyScale(rs->gf, B, rs->T, GfDiv(rs->gf, d, b), B);
            GfPolySub(rs->gf, T, rs->T, B, rs->T, C);
            //restore T2 to B
            memcpy(B, T2, rs->T * sizeof(*B));
            m++;
        }
    }

    *outSize = L + 1;
    memcpy(out, C, rs->T * sizeof(*out));

    if((L * 2) > rs->T)
        status = CFEC_UNCORRECTABLE;

leave:
    free(B);
    free(C);
    free(T);
    free(T2);
    return status;
#else

    DATA_TYPE *errLoc = calloc((rs->T + 1), sizeof(*errLoc));
    DATA_TYPE *newLoc = calloc((rs->T + 1), sizeof(*newLoc));
    DATA_TYPE *oldLoc = calloc((rs->T + 1), sizeof(*oldLoc));
    DATA_TYPE *tmpLoc = calloc((rs->T + 1), sizeof(*tmpLoc));

    if(!errLoc || !newLoc || !oldLoc || !tmpLoc)
    {
        status = CFEC_NO_RESOURCES;
        goto leave;
    }

    size_t newLocLen = 0;
    size_t errLocLen = 1;
    size_t oldLocLen = 1;

    errLoc[0] = 1;
    oldLoc[0] = 1;

    for(size_t i = 0; i < rs->T; i++)
    {
        DATA_TYPE delta = syndromes[i];
        for(size_t j = 1; j < errLocLen; j++)
        {
            delta = GfSub(rs->gf, delta, GfMul(rs->gf, errLoc[errLocLen - j - 1], syndromes[i - j]));
        }

        oldLoc[oldLocLen] = 0;
        oldLocLen++;

        if(delta != 0)
        {
            if(oldLocLen > errLocLen)
            {
                GfPolyScale(rs->gf, oldLoc, oldLocLen, delta, newLoc);
                newLocLen = oldLocLen;
                GfPolyScale(rs->gf, errLoc, errLocLen, GfInv(rs->gf, delta), oldLoc);
                oldLocLen = errLocLen;
                memcpy(errLoc, newLoc, newLocLen * sizeof(*errLoc));
                errLocLen = newLocLen;
            }
            //newLoc won't be used anymore, use it as a buffer
            GfPolyScale(rs->gf, oldLoc, oldLocLen, delta, newLoc);
            newLocLen = oldLocLen;
            memcpy(tmpLoc, errLoc, errLocLen * sizeof(*errLoc));
            errLocLen = GfPolyAdd(rs->gf, tmpLoc, errLocLen, newLoc, newLocLen, errLoc);
        }
    }

    size_t index = 0;
    for(size_t i = 0; i < errLocLen; i++)
    {
        if((index == 0) && (errLoc[i] == 0)) //drop leading zeros
            continue;

        out[index++] = errLoc[i];
    }

    if(((index - 1) << 1) > rs->T)
    {
        status = CFEC_UNCORRECTABLE;
        goto leave;
    }

    GfPolyInv(rs->gf, out, index);

    *outSize = index;

leave:
    free(errLoc);
    free(newLoc);
    free(oldLoc);
    free(tmpLoc);
    return status;
#endif
}

/**
 * @brief Calculates error magnitude (errata) polynomial and fix data
 * @param *rs RS instance
 * @param *data Input data block
 * @param *syn Syndrome polynomial
 * @param *locator Error locator polynomial
 * @param *position Error position vector
 * @param errCount Number of errors in position vector (also locator length - 1)
 * @return Status code
 */
static CFecStatus RsFix(struct CFecRS *rs, DATA_TYPE *data, 
    const DATA_TYPE *syndromes, const DATA_TYPE *locator, const DATA_TYPE *positions, size_t errCount)
{
    CFecStatus status = CFEC_OK;
    /*
     * This is based on Forney's algorithm.
     */

    //error evaluator polynomial Omega(x)
    DATA_TYPE *errorEvaluator = calloc((2 * rs->T + 2), sizeof(*errorEvaluator)); 
    //vector of error location values X_l, which are the zeros of the error locator polynomial Lambda(x) (passed in *locator)
    DATA_TYPE *errorLocation = malloc(rs->T * sizeof(*errorLocation));

    if(!errorEvaluator || !errorLocation)
    {
        status = CFEC_NO_RESOURCES;
        goto leave;
    }
    
    //find Omega(x)=S(x)*Lambda(x) mod x^(t+1), where Omega(x) is the error evaluator polynomial,
    //S(x) is the syndrome polynomial, Lambda(x) is the error locator polynomial,
    //and t is the number of errors (error positions)
    GfPolyMul(rs->gf, syndromes, rs->T, locator, errCount + 1, errorEvaluator);
    
    //S(x) and Lambda(x) are stored lowest-degree first, so the result would be lowest-degree first too.
    //Then, we need to invert this polynomial so the evaluation later works.
    //Also, we don't really want to do the polynomial division, but since we divide by x^(t+1), we just need to truncate
    //the coefficients. This is done by limiting the length passed to GfPolyInv().
    GfPolyInv(rs->gf, errorEvaluator, errCount + 1);

    //Construct error location vector X_l, which are the roots of the error locator polynomial Lambda(x)
    //The vector of error positions contains reversed indices, so we need to reverse them back
    //This means we find Lambda(X_l)=0, which is alpha^i. This is done in RsErrorPositions() to find i's for these zeros.
    for(size_t i = 0; i < errCount; i++)
    {
        errorLocation[i] = GfPowAlpha(rs->gf, rs->N - 1 - positions[i]);
    }

    /*
    Use Forney algorithm to find error values (magnitudes).
    It is typically defined as Y_l=-X_l^(1-c)*Omega(X_l^-1)/Lambda'(X_l^-1) when Omega(x) was calculated mod x^(t+1),
    but this requires the derivative calculation (Lambda'(x)).
    Different sources, e.g. Blahut, Algebraic Codes for Data Transmission, simplify this to omit the derivative:
    Y_l=X_l^(-c)*Omega(X_l^-1)/(prod for j!=l of (1-X_j*X_l^-1))
    The book does not include the X_l^(-c) term which is needed for FCR different than 0,
    but here it is already present: https://web.archive.org/web/20140630172526/http://web.stanford.edu/class/ee387/handouts/notes7.pdf
    Note that in this lecture Omega(x) was calculated mod x^2t, so this term is X_l^(1-c).
    Also, the product for j=l is equal to 0, so we skip this to not divide by zero.
    */
    for(size_t i = 0; i < errCount; i++)
    {
        //get X_l^(-1), which is used for further calculations
        DATA_TYPE errorLocationInv = GfInv(rs->gf, errorLocation[i]);

        //calculate the denominator Lambda'(x) and evaluate it at X_l^(-1)
        DATA_TYPE locatorPrime = 1; 
        for(size_t j = 0; j < errCount; j++)
        {
            if(j != i)
            {
                locatorPrime = GfMul(rs->gf, locatorPrime, GfSub(rs->gf, 1, GfMul(rs->gf, errorLocationInv, errorLocation[j])));
            }
        }

        //evaluate error evaluator Omega(x) at X_l^(-1)
        DATA_TYPE magnitude = GfPolyEval(rs->gf, errorEvaluator, errCount + 1, errorLocationInv);
        
        //multiply by X_i^(-c)=X_l^(c) to account for non-zero FCR
        //for FCR=0 this would not change the magnitude, so skip
        if(rs->fcr > 0)
            magnitude = GfMul(rs->gf, magnitude, GfPow(rs->gf, errorLocationInv, rs->fcr));

        //make sure the denominator is non-zero, so that we don't divide by zero
        if(0 == locatorPrime)
        {
            status = CFEC_UNCORRECTABLE;
            goto leave;
        }

        //divide by Lambda'(X_l^(-1)) to get the error magnitude
        magnitude = GfDiv(rs->gf, magnitude, locatorPrime);

        data[positions[i]] = GfSub(rs->gf, data[positions[i]], magnitude);
    }

leave:
    free(errorEvaluator);
    free(errorLocation);
    return status;
}


/**
 * @brief Check if syndromes indicate no error in the message
 * @param *syndromes Syndrome polynomial
 * @param size Syndrome polynomial size (buffer length)
 * @return True if all syndromes are zero and the message is correct
*/
static bool RsCheckSyndromes(const DATA_TYPE *syndromes, size_t size)
{
    for(size_t i = 0; i < size; i++) //check if all syndromes are 0, if so, the message is correct
    {
        if(syndromes[i] != 0)
        {
            return false;
        }
    }
    return true;
}

CFecStatus RsDecode(struct CFecRS *rs, DATA_TYPE *data, size_t size, size_t *fixed)
{
    if(!rs || !data)
        return CFEC_BAD_PARAMETER;

    CFecStatus status = CFEC_OK;
    bool complete = (size == (rs->N - rs->T)) ? true : false;

    if(size > (rs->N - rs->T))
        return CFEC_BAD_PARAMETER;
    
    DATA_TYPE *buffer = complete ? data : malloc((rs->N) * sizeof(*buffer));
    DATA_TYPE *syndromes = malloc((rs->T + 1) * sizeof(*syndromes));
    DATA_TYPE *locator = malloc((rs->T + 1) * sizeof(*locator));
    DATA_TYPE *positions = malloc((rs->T + 1) * sizeof(*positions));

    if(!syndromes || !locator || !positions || (!complete && !buffer))
    {
        status = CFEC_NO_RESOURCES;
        goto leave;
    }

    //if the message fully occupies the block, we can directly work on it
    //otherwise, we need to move the received parity symbols to the end of the block
    if(!complete)
    {
        memset(buffer, 0, rs->N * sizeof(*buffer));
        memcpy(buffer, data, size * sizeof(*buffer));
        memcpy(&buffer[rs->N - rs->T], &data[size], (rs->T) * sizeof(*buffer));
    }
    
    //calculate syndrome polynomial S(x)
    RsSyndromes(rs, buffer, rs->N, syndromes);

    //check if there is anything to correct
    if(RsCheckSyndromes(syndromes, rs->T))
        goto leave;

    //calculate error locator polynomial Lambda(x)
    size_t locatorSize = 0;
    status = RsErrorLocator(rs, syndromes, locator, &locatorSize); //calculate error locator polynomial

    if(CFEC_OK != status)
        goto leave;

    //calculate vector of error positions
    status = RsErrorPositions(rs, locator, locatorSize, positions);

    if(CFEC_OK != status)
        goto leave;

    status = RsFix(rs, buffer, syndromes, locator, positions, locatorSize - 1); 
    
    if(CFEC_OK != status)
        goto leave;

    RsSyndromes(rs, buffer, rs->N, syndromes); //calculate syndromes again to check if the message has been corrected successfully

    if(RsCheckSyndromes(syndromes, rs->T))
    {
        if(NULL != fixed)
            *fixed = locatorSize - 1;

        if(!complete)
            memcpy(data, buffer, size * sizeof(*buffer));
    }
    else
        status = CFEC_UNCORRECTABLE;

leave:
    free(syndromes);
    free(locator);
    free(positions);
    if(!complete)
        free(buffer);
    return status;
}

CFecStatus RsEncode(struct CFecRS *rs, DATA_TYPE *data, size_t size)
{
    if(!rs)
        return CFEC_BAD_PARAMETER;

    if(size > (rs->N - rs->T))
        return CFEC_BAD_PARAMETER;

    DATA_TYPE *buffer = calloc(rs->N, sizeof(*buffer));
    DATA_TYPE *divBuffer = malloc(rs->N * sizeof(*divBuffer));
    if(!buffer || !divBuffer)
    {
        free(buffer);
        free(divBuffer);
        return CFEC_NO_RESOURCES;
    }

    memcpy(buffer, data, size * sizeof(*buffer));
    //Divide message by the generator polynomial and copy the remainder polynomial (parity symbols) to immediately
    //succeed the message
    memcpy(&data[size], GfPolyDiv(rs->gf, buffer, rs->N, rs->generator, rs->T + 1, divBuffer), rs->T * sizeof(*data));
    //Since the systematic equation is C(x)=x^(n-k)M(x)-R(x), the remainder needs to be negated
    GfPolyNeg(rs->gf, &data[size], rs->T);

    free(buffer);
    free(divBuffer);
    return CFEC_OK;
}

static bool RsFindGfParams(size_t N, size_t *p, size_t *n)
{
    N = N + 1;
    *n = 0;
    bool found = false;

    if(0 == (N % 2))
    {
        *p = 2;
        found = true;
        while(0 == (N % 2))
        {
            (*n)++;
            N /= 2;
        }
    }

    for(size_t i = 3; (i * i) <= N; i += 2)
    {
        if(0 == (N % i))
        {
            if(found)
                return false;
            else
                found = true;
            *p = i;
            while(0 == (N % i))
            {
                N /= i;
                (*n)++;
            }
        }
    }

    if(N > 1)
    {
        if(found)
            return false;
        else
        {
            *p = N;
            *n = 1;
        }
    }

    return true;
}


CFecStatus RsInit(size_t N, size_t T, size_t fcr, size_t primitive, struct CFecRS **rsUser)
{
    struct CFecRS *rs = NULL;
    CFecStatus status = CFEC_OK;
    size_t p = 0, n = 0;

    if(!RsFindGfParams(N, &p, &n))
        return CFEC_BAD_PARAMETER;

    if(NULL == rsUser)
        return CFEC_BAD_PARAMETER;

    rs = malloc(sizeof(*rs) + (T + 1) * sizeof(*(rs->generator)));
    if(NULL == rs)
        return CFEC_NO_RESOURCES;

    status = GfInit(p, n, primitive, &(rs->gf));
    if(CFEC_OK != status)
    {
        free(rs);
        return status;
    }
    
    DATA_TYPE *buffer = malloc((T + 1) * sizeof(*buffer));
    if(NULL == buffer)
    {
        GfDeinit(rs->gf);
        free(rs);
        return CFEC_NO_RESOURCES;
    }
    memset(rs->generator, 0, (T + 1) * sizeof(*(rs->generator)));
    //Build the generator polynomial G(x)=prod 0<=i<T of (x-alpha^(i+fcr)) to be used for encoding
    rs->generator[0] = 1;
    for(size_t i = 0; i < T; i++)
    {
        memcpy(buffer, rs->generator, (i + 1) * sizeof(*buffer));
        GfPolyMul(rs->gf, buffer, i + 1, (DATA_TYPE[2]){1, GfNeg(rs->gf, GfPowAlpha(rs->gf, i + fcr))}, 2, rs->generator);
    }
    rs->T = T;
    rs->N = N;
    rs->fcr = fcr;
    *rsUser = rs;
    return CFEC_OK;
}

CFecStatus RsDeinit(struct CFecRS *rs)
{
    if(nullptr == rs)
        return CFEC_BAD_PARAMETER;

    if(nullptr != rs->gf)
        GfDeinit(rs->gf);

    free(rs);
    return CFEC_OK;
}