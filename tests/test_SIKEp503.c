/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: benchmarking/testing isogeny-based key encapsulation mechanism SIKEp503
*********************************************************************************************/ 

#include <stdio.h>
#include <string.h>
#include "test_extras.h"
#include "../src/P503/P503_api.h"


#define SCHEME_NAME    "SIKEp503"

#define Hy_crypto_kem_keypair         Hy_crypto_kem_keypair_SIKEp503
#define Hy_crypto_kem_enc             Hy_crypto_kem_enc_SIKEp503
#define Hy_crypto_kem_dec             Hy_crypto_kem_dec_SIKEp503
#include "test_sike.c"