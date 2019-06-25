/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: benchmarking/testing isogeny-based key encapsulation mechanism SIKEp434
*********************************************************************************************/ 

#include <stdio.h>
#include <string.h>
#include "test_extras.h"
#include "../src/P434/P434_api.h"


#define SCHEME_NAME    "SIKEp434"

#define Hy_crypto_kem_keypair         Hy_crypto_kem_keypair_SIKEp434
#define Hy_crypto_kem_enc             Hy_crypto_kem_enc_SIKEp434
#define Hy_crypto_kem_dec             Hy_crypto_kem_dec_SIKEp434

#include "test_sike.c"