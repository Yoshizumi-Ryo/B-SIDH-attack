

assert(IsPrime(p));
assert((p mod 4) eq 3);
assert(N_A gt N_B); 
assert(IsOdd(N_A));
assert(IsDivisibleBy(p+1,N_A) or IsDivisibleBy(p-1,N_A));
assert(IsDivisibleBy(p+1,N_B) or IsDivisibleBy(p-1,N_B));
assert((N_A-N_B)*N_B gt p); 


//the 8th primitive root of 1.====================
_<x>:=PolynomialRing(GF(p^2));
zeta_8:=RootsInSplittingField(x^4+1)[1][1];
//=======================================


//public construction.==================
//E_0: y^2=x^3-x=x(x-1)(x+1).
_<x>:=PolynomialRing(GF(p^2));
E_0:=EllipticCurve(x*(x-1)*(x+1));
assert(IsSupersingular(E_0));
E_0_4:=BaseChange(E_0,GF(p^4));

//take one basis of E_0[N_A]=(Z/N_A Z)^2.
P_A,Q_A:=ell_to_torsion_basis_2(E_0_4,N_A);
//take one basis of E_0[N_B]=(Z/N_B Z)^2.
P_B,Q_B:=ell_to_torsion_basis_2(E_0_4,N_B);
//========================================



//Bob calculates secretly.==================
coff_B:=Random(0,(N_B-1));
R_B:=E_0_4!(P_B+coff_B*Q_B);
assert (Order(R_B) eq N_B);
E_B,PA_EB,QA_EB:=elliptic_isogeny_1ptker(E_0_4,R_B,P_A,Q_A);

assert(Order(PA_EB) eq N_A);
assert(Order(QA_EB) eq N_A);

//Note that the following data are public.
//E_A,PB_EA,QB_EA,E_B,PA_EB,QA_EB;
//===========================================


//construction auxiliary poinsts.================
E_pr,alpha_P_A,alpha_Q_A:=construct_auxiliary_img_6(E_0_4,N_A,N_B,P_A,Q_A);
//===============================================



//main attack.================================
ker:=MainTorsionAttackKummer_2(E_0_4,E_B,E_pr,N_A,N_B,P_A,Q_A,PA_EB,QA_EB,alpha_P_A,alpha_Q_A,zeta_8);
//===============================================



//check if the attack succeed.=============
"";
"Result of the attack.",
(Order(ker) eq N_B) and 
WeilPairing(E_0_4!ker,E_0_4!R_B,N_B) eq 1;
//==========================================
