/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: elliptic curve and isogeny functions
*********************************************************************************************/


void xDBL(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t C24)
{ // Doubling of a Montgomery point in projective coordinates (X:Z).
  // Input: projective Montgomery x-coordinates P = (X1:Z1), where x1=X1/Z1 and Montgomery curve constants A+2C and 4C.
  // Output: projective Montgomery x-coordinates Q = 2*P = (X2:Z2).
    f2elm_t t0, t1;
    
    fp2sub(P->X, P->Z, t0);                         // t0 = X1-Z1
    fp2add(P->X, P->Z, t1);                         // t1 = X1+Z1
    fp2sqr_mont(t0, t0);                            // t0 = (X1-Z1)^2 
    fp2sqr_mont(t1, t1);                            // t1 = (X1+Z1)^2 
    fp2mul_mont(C24, t0, Q->Z);                     // Z2 = C24*(X1-Z1)^2   
    fp2mul_mont(t1, Q->Z, Q->X);                    // X2 = C24*(X1-Z1)^2*(X1+Z1)^2
    fp2sub(t1, t0, t1);                             // t1 = (X1+Z1)^2-(X1-Z1)^2 
    fp2mul_mont(A24plus, t1, t0);                   // t0 = A24plus*[(X1+Z1)^2-(X1-Z1)^2]
    fp2add(Q->Z, t0, Q->Z);                         // Z2 = A24plus*[(X1+Z1)^2-(X1-Z1)^2] + C24*(X1-Z1)^2
    fp2mul_mont(Q->Z, t1, Q->Z);                    // Z2 = [A24plus*[(X1+Z1)^2-(X1-Z1)^2] + C24*(X1-Z1)^2]*[(X1+Z1)^2-(X1-Z1)^2]
}


void xDBLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t C24, const int e)
{ // Computes [2^e](X:Z) on Montgomery curve with projective constant via e repeated doublings.
  // Input: projective Montgomery x-coordinates P = (XP:ZP), such that xP=XP/ZP and Montgomery curve constants A+2C and 4C.
  // Output: projective Montgomery x-coordinates Q <- (2^e)*P.
    int i;
    
    copy_words((digit_t*)P, (digit_t*)Q, 2*2*NWORDS_FIELD);

    for (i = 0; i < e; i++) {
        xDBL(Q, Q, A24plus, C24);
    }
}

#if (OALICE_BITS % 2 == 1)

void get_2_isog(const point_proj_t P, f2elm_t A, f2elm_t C)
{ // Computes the corresponding 2-isogeny of a projective Montgomery point (X2:Z2) of order 2.
  // Input:  projective point of order two P = (X2:Z2).
  // Output: the 2-isogenous Montgomery curve with projective coefficients A/C.
    
    fp2sqr_mont(P->X, A);                           // A = X2^2
    fp2sqr_mont(P->Z, C);                           // C = Z2^2
    fp2sub(C, A, A);                                // A = Z2^2 - X2^2
}


void eval_2_isog(point_proj_t P, point_proj_t Q)
{ // Evaluates the isogeny at the point (X:Z) in the domain of the isogeny, given a 2-isogeny phi.
  // Inputs: the projective point P = (X:Z) and the 2-isogeny kernel projetive point Q = (X2:Z2).
  // Output: the projective point P = phi(P) = (X:Z) in the codomain. 
    f2elm_t t0, t1, t2, t3;
    
    fp2add(Q->X, Q->Z, t0);                         // t0 = X2+Z2
    fp2sub(Q->X, Q->Z, t1);                         // t1 = X2-Z2
    fp2add(P->X, P->Z, t2);                         // t2 = X+Z
    fp2sub(P->X, P->Z, t3);                         // t3 = X-Z
    fp2mul_mont(t0, t3, t0);                        // t0 = (X2+Z2)*(X-Z)
    fp2mul_mont(t1, t2, t1);                        // t1 = (X2-Z2)*(X+Z)
    fp2add(t0, t1, t2);                             // t2 = (X2+Z2)*(X-Z) + (X2-Z2)*(X+Z)
    fp2sub(t0, t1, t3);                             // t3 = (X2+Z2)*(X-Z) - (X2-Z2)*(X+Z)
    fp2mul_mont(P->X, t2, P->X);                    // Xfinal
    fp2mul_mont(P->Z, t3, P->Z);                    // Zfinal
}

#endif



void get_4_isog_ED_jinv(const point_proj_t P, f2elm_t C, f2elm_t D)
{ // Computes the corresponding 4-isogeny of a projective Montgomery point (X4:Z4) of order 4.
  // Input:  projective point of order four P = (X4:Z4).
  // Output: the 4-isogenous Edwards curve with projective coefficients C=C+D/D=C-D for computing j-invariant 

	f2elm_t t1;

	fp2sqr_mont(P->X, C); 							// C = (1/4)(Y4+Z4)^2
	fp2sqr_mont(P->Z, D); 							// t1 = (1/4)(Y4-Z4)^2
	fp2sqr_mont(C, C);							    // C = (1/16)(Y4+Z4)^4 // C
	fp2sqr_mont(D, D);							    // D = (1/16)(Y4-Z4)^4 // C-D
	fp2sub(C, D, t1);								// t1 = D
	fp2add(C, t1, C);								// C = C+D
}



void get_4_isog_ED(const point_proj_t P, f2elm_t C, f2elm_t D, f2elm_t coeff)
{ // Computes the corresponding 4-isogeny of a projective Edwards point (X4:Z4) of order 4.
  // Input:  projective point of order four P = (X4:Z4).
  // Output: the 4-isogenous Montgomery curve with projective coefficients C'=A+2C D'=4C and the coefficient
  //         that are used to evaluate the isogeny at a point in eval_4_isog().

	fp2sqr_mont(P->X, C); 							// C = 1/4(Y4+Z4)^2
	fp2sqr_mont(P->Z, D);							// t0 = 1/4(Y4-Z4)^2
	fp2add(C, C, coeff); 						    // coeff = 1/2(Y4+Z4)^2
	fp2add(coeff, coeff, coeff); 				    // coeff = (Y4+Z4)^2
	fp2sqr_mont(C, C); 								// C = 1/4(Y4+Z4)^4 // C
	fp2sqr_mont(D, D);								// t0 = 1/4(Y4-Z4)^4 //C-D

}


void eval_4_isog_ED(point_proj_ED_t P, f2elm_t coeff, point_proj_t M)
{ // Evaluates the isogeny at the point (X:Z) in the domain of the isogeny, given a 4-isogeny phi defined 
  // by the coefficient in coeff (computed in the function get_4_isog()).
  // Inputs: the coefficient defining the isogeny, kernel point P on an Ewards curve, and the projective point M = (X:Z).
  // Output: the projective point M = phi(M) = (X:Z) in the codomain. 
	f2elm_t t0, t1, t2, t3;

	fp2add(M->X, M->Z, t1); 							// Z
	fp2sub(M->X, M->Z, t0);								// Y
	fp2mul_mont(t0, t1, t2); 							// t2 = YZ
	fp2mul_mont(t2, coeff, t2); 						// t2 = YZ*(Y4+Z4)^2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             yuhtjj j Z*(Y4+Z4)^2
	fp2mul_mont(t0, P->Z, t0); 							// t0 = Y*Z4
	fp2mul_mont(t1, P->Y, t1); 							// t1 = Z*Y4
	fp2add(t0, t1, t3); 								// t3 = (YZ4+ZY4)
	fp2sub(t0, t1, t0); 								// t0 = (YZ4-ZY4)
	fp2sqr_mont(t3, t3); 								// t3 = (YZ4+ZY4)^2
	fp2sqr_mont(t0, t0); 								// t0 = (YZ4-ZY4)^2	
    fp2add(t2, t0, M->X);								// M->X = YZ*(Y4+Z4)^2 +  (YZ4-ZY4)^2
	fp2mul_mont(M->X, t3, M->X);						// M->X = (YZ4+ZY4)^2*(YZ*(Y4+Z4)^2 +  (YZ4-ZY4)^2)
	fp2sub(t3, t2, M->Z);
	fp2mul_mont(M->Z, t0, M->Z);

}



void xTPL_ED(const point_proj_t P, point_proj_t Q, const f2elm_t C, const f2elm_t D)
{ // Tripling of a Montgomery point in projective coordinates (X:Z) using Edwards' tripling formula.
  // Input: projective Montgomery x-coordinates P = (X:Z), where x=X/Z and Edwards curve constants C and D.
  // Output: projective Montgomery x-coordinates Q = 3*P = (X3:Z3).
	f2elm_t t0, t1, t2, t3, t4, t5;

	fp2sub(P->X, P->Z, t0);
	fp2add(P->X, P->Z, t1);    
    fp2sub(t0, t1, t2); 						    // t2 = Y-Z
	fp2sqr_mont(t0, t0); 							// t0 = Y^2
	fp2sqr_mont(t1, t1); 							// t1 = Z^2
	fp2sqr_mont(t2, t3);							// t3 = (Y-Z)^2
	fp2sub(t3, t0, t3); 							// t3 = Z^2-2YZ
	fp2sub(t3, t1, t3); 							// t3 = -2YZ
	fp2mul_mont(t0, D, t4); 						// t4 = DY^2
	fp2mul_mont(t1, C, t5);							// t5 = CZ^2
	fp2sub(t4, t5, t2);								// t6 = DY^2-CZ^2
	fp2mul_mont(t2, t3, t2);						// t6 = -2YZ(DY^2-CZ^2)
	fp2mul_mont(t4, t0, t4);						// t4 = DY^4
	fp2mul_mont(t5, t1, t5);						// t5 = CZ^4
	fp2sub(t4, t5, t4);								// t4 = DY^4-CZ^4
	fp2add(t2, t4, t0);								// t0 = DY^2(Y^2+2ZY)-CZ^2(-2YZ+Z^2)
	fp2sub(t4, t2, t1);								// t1 = DY^2(Y^2-2ZY)+CZ^2(2YZ-Z^2)
	fp2sqr_mont(t0, t0);
	fp2sqr_mont(t1, t1);
	fp2mul_mont(t0, P->X, Q->X);					// F = (DY^2(Y^2+2ZY)-CZ^2(-2YZ+Z^2))^2(Y+Z)
	fp2mul_mont(t1, P->Z, Q->Z);					// G = (DY^2(Y^2-2ZY)+CZ^2(2YZ+Z^2))^2(Y-Z)

}




void xTPLe_ED(const point_proj_t P, point_proj_t Q, const f2elm_t C, const f2elm_t D, const int e)
{ // Computes [3^e](X:Z) on Montgomery curve with projective constant via e repeated triplings.
  // Input: projective Montgomery x-coordinates P = (XP:ZP), such that xP=XP/ZP and Montgomery curve constants A24plus = A+2C and A24minus = A-2C.
  // Output: projective Montgomery x-coordinates Q <- (3^e)*P.
	int i;

	copy_words((digit_t*)P, (digit_t*)Q, 2 * 2 * NWORDS_FIELD);

	for (i = 0; i < e; i++) {
		xTPL_ED(Q, Q, C, D);
	}
}




void get_3_isog_ED(const point_proj_ED_t P, f2elm_t C, f2elm_t D)
{ // Computes the corresponding 3-isogeny of a projective Edwards point (Y3:Z3) of order 3.
  // Input:  projective point of order three P = (Y3:Z3).
  // Output: the 3-isogenous Edwards curve with projective coefficient D/C. 
	f2elm_t t0, t1, t2, t3, t4;

	fp2add(P->Y, P->Z, t0);							// t0 = Y3+Z3
	fp2sqr_mont(t0, t0);							// t0 = (Y3+Z3)^2 = Y3^2+2Y3Z3+Z3^2
	fp2sqr_mont(P->Y, t1);							// t1 = Y3^2
	fp2sqr_mont(P->Z, t2);							// t2 = Z3^2
	fp2sub(t0, t1, t3);								// t3 = Z3^2+2Y3Z3
	fp2sub(t0, t2, t4);								// t4 = Y3^2+2Y3Z3
	fp2add(t0, t3, C); 								// C = Y3^2+4Y3Z3+2Z3^2			
	fp2add(t2, t2, t2); 							// t2 = 2Z3^2	
	fp2add(t2, C, C);								// C = Y3^2+4Y3Z3+4Z3^2 = (Y3+2Z3)^2
	fp2add(t0, t4, D); 								// D = Z3^2+4Y3Z3+2Y3^2	
	fp2add(t1, t1, t1);								// t1 = 2Y3^2
	fp2add(D, t1, D); 								// D = Z3^2+4Y3Z3+4Y3^2 = (Z3+2Y3)^2
	fp2mul_mont(C, t4, C);							// C = (Y3+2Z3)^2*(Y3^2+2Y3Z3)	
	fp2mul_mont(D, t3, D);	
}



void eval_3_isog_ED(point_proj_ED_t P,  point_proj_t M)
{ // Computes the 3-isogeny R=phi(X:Z), given projective point (Y3:Z3) of order 3 on an Edwards curve and 
  // a point M 
  // Inputs: projective points P = (Y3:Z3) and M = (X:Z).
  // Output: the projective point M <- phi(M) = (X:Z). 
	f2elm_t t0, t1, t2;

    fp2sub(M->X, M->Z, t1);                         // points to edwards 
	fp2add(M->X, M->Z, t0);
	fp2mul_mont(t0, P->Y, t0); 					    // t0 = ZY3
	fp2mul_mont(t1, P->Z, t1); 					    // t1 = YZ3
	fp2add(t0, t1, t2);								// t2 = ZY3+YZ3
	fp2sqr_mont(t2, t2); 							// t2 = (ZY3+YZ3)^2
	fp2mul_mont(t2, M->X, M->X); 				    // t2 = (ZY3+YZ3)^2*(Y+Z)
	fp2sub(t0, t1, t0);
	fp2sqr_mont(t0, t0); 							// t0 = (ZY3-YZ3)^2
	fp2mul_mont(t0, M->Z, M->Z); 				    // t0 = (ZY3-YZ3)^2*(Z-Y)

}


void inv_3_way(f2elm_t z1, f2elm_t z2, f2elm_t z3)
{ // 3-way simultaneous inversion
  // Input:  z1,z2,z3
  // Output: 1/z1,1/z2,1/z3 (override inputs).
    f2elm_t t0, t1, t2, t3;

    fp2mul_mont(z1, z2, t0);                      // t0 = z1*z2
    fp2mul_mont(z3, t0, t1);                      // t1 = z1*z2*z3
    fp2inv_mont(t1);                              // t1 = 1/(z1*z2*z3)
    fp2mul_mont(z3, t1, t2);                      // t2 = 1/(z1*z2) 
    fp2mul_mont(t2, z2, t3);                      // t3 = 1/z1
    fp2mul_mont(t2, z1, z2);                      // z2 = 1/z2
    fp2mul_mont(t0, t1, z3);                      // z3 = 1/z3
    fp2copy(t3, z1);                              // z1 = 1/z1
}


void get_A(const f2elm_t xP, const f2elm_t xQ, const f2elm_t xR, f2elm_t A)
{ // Given the x-coordinates of P, Q, and R, returns the value A corresponding to the Montgomery curve E_A: y^2=x^3+A*x^2+x such that R=Q-P on E_A.
  // Input:  the x-coordinates xP, xQ, and xR of the points P, Q and R.
  // Output: the coefficient A corresponding to the curve E_A: y^2=x^3+A*x^2+x.
    f2elm_t t0, t1, one = {0};
    
    fpcopy((digit_t*)&Montgomery_one, one[0]);
    fp2add(xP, xQ, t1);                           // t1 = xP+xQ
    fp2mul_mont(xP, xQ, t0);                      // t0 = xP*xQ
    fp2mul_mont(xR, t1, A);                       // A = xR*t1
    fp2add(t0, A, A);                             // A = A+t0
    fp2mul_mont(t0, xR, t0);                      // t0 = t0*xR
    fp2sub(A, one, A);                            // A = A-1
    fp2add(t0, t0, t0);                           // t0 = t0+t0
    fp2add(t1, xR, t1);                           // t1 = t1+xR
    fp2add(t0, t0, t0);                           // t0 = t0+t0
    fp2sqr_mont(A, A);                            // A = A^2
    fp2inv_mont(t0);                              // t0 = 1/t0
    fp2mul_mont(A, t0, A);                        // A = A*t0
    fp2sub(A, t1, A);                             // Afinal = A-t1
}




void j_inv_ED(const f2elm_t C, const f2elm_t D, f2elm_t jinv)
{ // Computes the j-invariant of a Edwards curve with projective constant.
  // Input: C+D, C-D in GF(p^2).
  // Output: j-invariant of the Edwards curve x^2+y^2=1+(D/C)x^2y^2
	f2elm_t t0, t1, t2, t3, t4, t5;

	fp2sqr_mont(C, t0); 							// t0 = (C+D)^2
	fp2sqr_mont(D, t1); 							// t1 = (C-D)^2
	fp2sub(t0, t1, t0); 							// t0 = 4CD
	fp2add(t0, t0, t2);								// t2 = 8CD
	fp2add(t2, t2, t2);								// t2 = 16CD
	fp2add(t1, t2, t2); 							// t2 = (C-D)^2+16CD
	fp2add(t2, t2, t2); 							// t2 = 2(C-D)^2+32CD
	fp2add(t2, t2, t2); 							// t2 = 4(C-D)^2+64CD
	fp2sqr_mont(t2, jinv);							// jinv = (4(C-D)^2+64CD)^2					
	fp2mul_mont(t2, jinv, jinv);					// jinv = (4(C-D)^2+64CD)^3
	fp2sqr_mont(t1, t1);							// t1 = (C-D)^4
	fp2mul_mont(t0, t1, t1);						// t1 = 4CD(C-D)^4
	fp2inv_mont(t1); 								// t1 = 1/4CD(C-D)^4
	fp2mul_mont(t1, jinv, jinv);
}


void xDBLADD(point_proj_t P, point_proj_t Q, const f2elm_t xPQ, const f2elm_t A24)
{ // Simultaneous doubling and differential addition.
  // Input: projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, affine difference xPQ=x(P-Q) and Montgomery curve constant A24=(A+2)/4.
  // Output: projective Montgomery points P <- 2*P = (X2P:Z2P) such that x(2P)=X2P/Z2P, and Q <- P+Q = (XQP:ZQP) such that = x(Q+P)=XQP/ZQP. 
    f2elm_t t0, t1, t2;

    fp2add(P->X, P->Z, t0);                         // t0 = XP+ZP
    fp2sub(P->X, P->Z, t1);                         // t1 = XP-ZP
    fp2sqr_mont(t0, P->X);                          // XP = (XP+ZP)^2
    fp2sub(Q->X, Q->Z, t2);                         // t2 = XQ-ZQ
    fp2correction(t2);
    fp2add(Q->X, Q->Z, Q->X);                       // XQ = XQ+ZQ
    fp2mul_mont(t0, t2, t0);                        // t0 = (XP+ZP)*(XQ-ZQ)
    fp2sqr_mont(t1, P->Z);                          // ZP = (XP-ZP)^2
    fp2mul_mont(t1, Q->X, t1);                      // t1 = (XP-ZP)*(XQ+ZQ)
    fp2sub(P->X, P->Z, t2);                         // t2 = (XP+ZP)^2-(XP-ZP)^2
    fp2mul_mont(P->X, P->Z, P->X);                  // XP = (XP+ZP)^2*(XP-ZP)^2
    fp2mul_mont(t2, A24, Q->X);                     // XQ = A24*[(XP+ZP)^2-(XP-ZP)^2]
    fp2sub(t0, t1, Q->Z);                           // ZQ = (XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)
    fp2add(Q->X, P->Z, P->Z);                       // ZP = A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2
    fp2add(t0, t1, Q->X);                           // XQ = (XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)
    fp2mul_mont(P->Z, t2, P->Z);                    // ZP = [A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2]*[(XP+ZP)^2-(XP-ZP)^2]
    fp2sqr_mont(Q->Z, Q->Z);                        // ZQ = [(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
    fp2sqr_mont(Q->X, Q->X);                        // XQ = [(XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)]^2
    fp2mul_mont(Q->Z, xPQ, Q->Z);                   // ZQ = xPQ*[(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
}


static void swap_points(point_proj_t P, point_proj_t Q, const digit_t option)
{ // Swap points.
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
    digit_t temp;
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++) {
        temp = option & (P->X[0][i] ^ Q->X[0][i]);
        P->X[0][i] = temp ^ P->X[0][i]; 
        Q->X[0][i] = temp ^ Q->X[0][i]; 
        temp = option & (P->Z[0][i] ^ Q->Z[0][i]);
        P->Z[0][i] = temp ^ P->Z[0][i]; 
        Q->Z[0][i] = temp ^ Q->Z[0][i]; 
        temp = option & (P->X[1][i] ^ Q->X[1][i]);
        P->X[1][i] = temp ^ P->X[1][i]; 
        Q->X[1][i] = temp ^ Q->X[1][i]; 
        temp = option & (P->Z[1][i] ^ Q->Z[1][i]);
        P->Z[1][i] = temp ^ P->Z[1][i]; 
        Q->Z[1][i] = temp ^ Q->Z[1][i]; 
    }
}


static void LADDER3PT(const f2elm_t xP, const f2elm_t xQ, const f2elm_t xPQ, const digit_t* m, const unsigned int AliceOrBob, point_proj_t R, const f2elm_t A)
{
    point_proj_t R0 = {0}, R2 = {0};
    f2elm_t A24 = {0};
    digit_t mask;
    int i, nbits, bit, swap, prevbit = 0;

    if (AliceOrBob == ALICE) {
        nbits = OALICE_BITS;
    } else {
        nbits = OBOB_BITS - 1;
    }

    // Initializing constant
    fpcopy((digit_t*)&Montgomery_one, A24[0]);
    fp2add(A24, A24, A24);
    fp2add(A, A24, A24);
    fp2div2(A24, A24);  
    fp2div2(A24, A24);  // A24 = (A+2)/4

    // Initializing points
    fp2copy(xQ, R0->X);
    fpcopy((digit_t*)&Montgomery_one, (digit_t*)R0->Z);
    fp2copy(xPQ, R2->X);
    fpcopy((digit_t*)&Montgomery_one, (digit_t*)R2->Z);
    fp2copy(xP, R->X);
    fpcopy((digit_t*)&Montgomery_one, (digit_t*)R->Z);
    fpzero((digit_t*)(R->Z)[1]);

    // Main loop
    for (i = 0; i < nbits; i++) {
        bit = (m[i >> LOG2RADIX] >> (i & (RADIX-1))) & 1;
        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (digit_t)swap;

        swap_points(R, R2, mask);
        xDBLADD(R0, R2, R->X, A24);
        fp2mul_mont(R2->X, R->Z, R2->X);
    }
    swap = 0 ^ prevbit;
    mask = 0 - (digit_t)swap;
    swap_points(R, R2, mask);
}