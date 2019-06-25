/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: benchmarking/testing isogeny-based key exchange SIDHp434
*********************************************************************************************/ 

#include <stdio.h>
#include <string.h>
#include "test_extras.h"
#include "../src/P434/P434_api.h"


#define SCHEME_NAME    "SIDHp434"

#define random_mod_order_A            random_mod_order_A_SIDHp434
#define random_mod_order_B            random_mod_order_B_SIDHp434


#define Hy_KeyGeneration_A	          Hy_KeyGeneration_A_SIDHp434
#define Hy_KeyGeneration_B	          Hy_KeyGeneration_B_SIDHp434
#define Hy_SecretAgreement_A		  Hy_SecretAgreement_A_SIDHp434
#define Hy_SecretAgreement_B		  Hy_SecretAgreement_B_SIDHp434

#include "test_sidh.c"