/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: supersingular isogeny parameters and generation of functions for P503
*********************************************************************************************/  

#include "P503_api.h" 
#include "P503_internal.h"


// Encoding of field elements, elements over Z_order, elements over GF(p^2) and elliptic curve points:
// --------------------------------------------------------------------------------------------------
// Elements over GF(p) and Z_order are encoded with the least significant octet (and digit) located at the leftmost position (i.e., little endian format). 
// Elements (a+b*i) over GF(p^2), where a and b are defined over GF(p), are encoded as {a, b}, with a in the least significant position.
// Elliptic curve points P = (x,y) are encoded as {x, y}, with x in the least significant position. 
// Internally, the number of digits used to represent all these elements is obtained by approximating the number of bits to the immediately greater multiple of 32.
// For example, a 503-bit field element is represented with Ceil(503 / 64) = 8 64-bit digits or Ceil(503 / 32) = 16 32-bit digits.

//
// Curve isogeny system "SIDHp503". Base curve: Montgomery curve By^2 = Cx^3 + Ax^2 + Cx defined over GF(p503^2), where A=6, B=1, C=1 and p503 = 2^250*3^159-1
//
         
const uint64_t p503[NWORDS64_FIELD]              = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xABFFFFFFFFFFFFFF, 
                                                     0x13085BDA2211E7A0, 0x1B9BF6C87B7E7DAF, 0x6045C6BDDA77A4D0, 0x004066F541811E1E };
const uint64_t p503p1[NWORDS64_FIELD]            = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0xAC00000000000000,
                                                     0x13085BDA2211E7A0, 0x1B9BF6C87B7E7DAF, 0x6045C6BDDA77A4D0, 0x004066F541811E1E };
const uint64_t p503x2[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x57FFFFFFFFFFFFFF,
                                                     0x2610B7B44423CF41, 0x3737ED90F6FCFB5E, 0xC08B8D7BB4EF49A0, 0x0080CDEA83023C3C }; 
// Order of Alice's subgroup
const uint64_t Alice_order[NWORDS64_ORDER]       = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0400000000000000 }; 
// Order of Bob's subgroup
const uint64_t Bob_order[NWORDS64_ORDER]         = { 0xC216F6888479E82B, 0xE6FDB21EDF9F6BC4, 0x1171AF769DE93406, 0x1019BD5060478798 };
// Alice's generator values {XPA0 + XPA1*i, XQA0 + XQA1*i, XRA0 + XRA1*i} in GF(p503^2), expressed in Montgomery representation
const uint64_t A_gen[6*NWORDS64_FIELD]           = { 0x5D083011589AD893, 0xADFD8D2CB67D0637, 0x330C9AC34FFB6361, 0xF0D47489A2E805A2,
                                                     0x27E2789259C6B8DC, 0x63866A2C121931B9, 0x8D4C65A7137DCF44, 0x003A183AE5967B3F,   // XPA0
                                                     0x7E3541B8C96D1519, 0xD3ADAEEC0D61A26C, 0xC0A2219CE7703DD9, 0xFF3E46658FCDBC52,
                                                     0xD5B38DEAE6E196FF, 0x1AAC826364956D58, 0xEC9F4875B9A5F27A, 0x001B0B475AB99843,   // XPA1
                                                     0x4D83695107D03BAD, 0x221F3299005E2FCF, 0x78E6AE22F30DECF2, 0x6D982DB5111253E4,
                                                     0x504C80A8AB4526A8, 0xEFD0C3AA210BB024, 0xCB77483501DC6FCF, 0x001052544A96BDF3,   // XQA0
                                                     0x0D74FE3402BCAE47, 0xDF5B8CDA832D8AED, 0xB86BCF06E4BD837E, 0x892A2933A0FA1F63,
                                                     0x9F88FC67B6CCB461, 0x822926EA9DDA3AC8, 0xEAC8DDE5855425ED, 0x000618FE6DA37A80,   // XQA1
                                                     0x6B6F4A4F786CF310, 0xB019D444BDCFDBE3, 0xA14CB06680607834, 0xCB0D5582E7E6E60A, 
                                                     0xBA4EE8771667E241, 0xE42A114FCB12E5FF, 0x9A0C074E275BCD98, 0x001871329B28689E,   // XRA0
                                                     0x685FA8378513FE76, 0x84E8FC1785E8BDF0, 0x8A380F177CB7C1B7, 0x2227464F4F812C94, 
                                                     0x117C94A81A90C279, 0x193D73132FB9FA28, 0x99335336F192C9EB, 0x0006C2FE778A34FD }; // XRA1
// Bob's generator values {XPB0, XQB0, XRB0 + XRB1*i} in GF(p503^2), expressed in Montgomery representation
const uint64_t B_gen[6*NWORDS64_FIELD]           = { 0xDF630FC5FB2468DB, 0xC30C5541C102040E, 0x3CDC9987B76511FC, 0xF54B5A09353D0CDD, 
                                                     0x3ADBA8E00703C42F, 0x8253F9303DDC95D0, 0x62D30778763ABFD7, 0x001CD00FB581CD55,   // XPB0
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,   // XPB1
                                                     0x2E3457A12B429261, 0x311F94E89627DCF8, 0x5B71C98FD1DB73F6, 0x3671DB7DCFC21541, 
                                                     0xB6D1484C9FE0CF4F, 0x19CD110717356E35, 0xF4F9FB00AC9919DF, 0x0035BC124D38A70B,   // XQB0
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,   // XQB1
                                                     0x2E08BB99413D2952, 0xD3021467CD088D72, 0x21017AF859752245, 0x26314ED8FFD9DE5C, 
                                                     0x4AF43C73344B6686, 0xCFA1F91149DF0993, 0xF327A95365587A89, 0x000DBF54E03D3906,   // XRB0
                                                     0xFC1FC00CBD0A0CFB, 0x66C29FB284B491A9, 0x7FBED0B26D7F18E0, 0x9C02361066F4C67D, 
                                                     0x2D83758DD0B19E6F, 0x0827029DB2CDC1D7, 0x58700A85FAD38A71, 0x0006FA256B614AD2 }; // XRB1
// Montgomery constant Montgomery_R2 = (2^512)^2 mod p503
const uint64_t Montgomery_R2[NWORDS64_FIELD]     = { 0x5289A0CF641D011F, 0x9B88257189FED2B9, 0xA3B365D58DC8F17A, 0x5BC57AB6EFF168EC,
                                                     0x9E51998BD84D4423, 0xBF8999CBAC3B5695, 0x46E9127BCE14CDB6, 0x003F6CFCE8B81771 };                                                   
// Value one in Montgomery representation 
const uint64_t Montgomery_one[NWORDS64_FIELD]    = { 0x00000000000003F9, 0x0000000000000000, 0x0000000000000000, 0xB400000000000000, 
                                                     0x63CB1A6EA6DED2B4, 0x51689D8D667EB37D, 0x8ACD77C71AB24142, 0x0026FBAEC60F5953 };


// Fixed parameters for isogeny tree computation
const unsigned int strat_Alice[MAX_Alice-1] = { 
61, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 
4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 
1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 29, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 
1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 13, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 
1, 1, 2, 1, 1, 5, 4, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1 };

const unsigned int strat_Bob[MAX_Bob-1] = { 
71, 38, 21, 13, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 5, 4, 2, 1, 1, 2, 1, 
1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 17, 9, 
5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 
1, 4, 2, 1, 1, 2, 1, 1, 33, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 
2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 1, 2, 
1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1 };
           
// Setting up macro defines and including GF(p), GF(p^2), curve, isogeny and kex functions
#define fpcopy                        fpcopy503
#define fpzero                        fpzero503
#define fpadd                         fpadd503
#define fpsub                         fpsub503
#define fpneg                         fpneg503
#define fpdiv2                        fpdiv2_503
#define fpcorrection                  fpcorrection503
#define fpmul_mont                    fpmul503_mont
#define fpsqr_mont                    fpsqr503_mont
#define fpinv_mont                    fpinv503_mont
#define fpinv_chain_mont              fpinv503_chain_mont
#define fpinv_mont_bingcd             fpinv503_mont_bingcd
#define fp2copy                       fp2copy503
#define fp2zero                       fp2zero503
#define fp2add                        fp2add503
#define fp2sub                        fp2sub503
#define fp2neg                        fp2neg503
#define fp2div2                       fp2div2_503
#define fp2correction                 fp2correction503
#define fp2mul_mont                   fp2mul503_mont
#define fp2sqr_mont                   fp2sqr503_mont
#define fp2inv_mont                   fp2inv503_mont
#define fp2inv_mont_bingcd            fp2inv503_mont_bingcd
#define fpequal_non_constant_time     fpequal503_non_constant_time
#define mp_add_asm                    mp_add503_asm
#define mp_subx2_asm                  mp_sub503x2_asm
#define mp_dblsubx2_asm               mp_dblsub503x2_asm
#define Hy_crypto_kem_keypair         Hy_crypto_kem_keypair_SIKEp503
#define Hy_crypto_kem_enc             Hy_crypto_kem_enc_SIKEp503
#define Hy_crypto_kem_dec             Hy_crypto_kem_dec_SIKEp503
#define random_mod_order_A            random_mod_order_A_SIDHp503
#define random_mod_order_B            random_mod_order_B_SIDHp503
#define Hy_KeyGeneration_A	          Hy_KeyGeneration_A_SIDHp503
#define Hy_KeyGeneration_B		      Hy_KeyGeneration_B_SIDHp503
#define Hy_SecretAgreement_A		  Hy_SecretAgreement_A_SIDHp503
#define Hy_SecretAgreement_B		  Hy_SecretAgreement_B_SIDHp503
#include "../fpx.c"
#include "../ec_isogeny.c"
#include "../sidh.c"    
#include "../sike.c"