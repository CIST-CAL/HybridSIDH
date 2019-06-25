/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: benchmarking/testing isogeny-based key exchange SIDHp610
*********************************************************************************************/ 

#include <stdio.h>
#include <string.h>
#include "test_extras.h"
#include "../src/P610/P610_api.h"


#define SCHEME_NAME    "SIDHp610"

#define random_mod_order_A            random_mod_order_A_SIDHp610
#define random_mod_order_B            random_mod_order_B_SIDHp610


#define Hy_KeyGeneration_A	          Hy_KeyGeneration_A_SIDHp610
#define Hy_KeyGeneration_B	          Hy_KeyGeneration_B_SIDHp610
#define Hy_SecretAgreement_A		  Hy_SecretAgreement_A_SIDHp610
#define Hy_SecretAgreement_B		  Hy_SecretAgreement_B_SIDHp610
#include "test_sidh.c"