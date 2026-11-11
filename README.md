# CFEC
Universal Forward Error Correction library written in pure C.

This library is in an early stage of development. Started this for fun and to implement some non-typical FEC algorithms.

For now, the library supports:
* Reed-Solomon in GF(2^k), GF(p), and GF(p^n) - to some extent - the appropriate field is chosen based on the provided N (block size).
