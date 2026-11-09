#ifndef CFEC_H_
#define CFEC_H_

/**
 * @brief CFEC return status codes
 */
typedef enum CFecStatus
{
    CFEC_OK = 0, /**< Operation successful */
    CFEC_NO_RESOURCES = 1, /**< Out of resources (e.g., memory) */
    CFEC_BAD_PARAMETER = 2, /**< Bad parameter */
    CFEC_UNCORRECTABLE = 3, /**< FEC library reports uncorrectable errors in data */
} CFecStatus;

#endif