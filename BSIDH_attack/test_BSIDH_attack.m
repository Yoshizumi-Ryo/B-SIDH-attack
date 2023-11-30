//これはあまり書き換えるな!!

//WLOG, we assume N_A>N_B.
//Over F_(p^4), we use Alice's model of elliptic curve in P^2.

//we consider theta structure given by Thomae formula for Legendre form y^2=x(x-1)(x-lmd).
//---------------


K:=GF(p);
_<t>:=PolynomialRing(GF(p^4));

/*parametor setting.===========================

We assume N_A>N_B.
Then, 
N_A: Degree of isogney using attack.
N_B: Order of torsion.
we assume N_A is odd.*/


/*
N_A:=3^3*5; 
N_B:=2^4;   
*/

assert(N_A gt N_B); 
assert(IsOdd(N_A));


Prime_Fac_N_A:=Seqset(fatoriztion_seq(N_A));
precomp_for_N_A:=AssociativeArray();

for l in Prime_Fac_N_A do
  precomp_for_N_A[l]:=AssociativeArray();
  Mat_F:=const_Mat_F(l);
  r,index_t,index_j:=const_index_t_j_2(l,Mat_F); //need.time
  precomp_for_N_A[l]["Mat_F"]:=Mat_F;
  precomp_for_N_A[l]["r"]:=r;
  precomp_for_N_A[l]["index_t"]:=index_t;
  precomp_for_N_A[l]["index_j"]:=index_j;
 end for;


//================================


//public construction.==================

//E_0: y^2=x^3-x=x(x-1)(x+1).
lmd_0:=K!(-1);
_,lv22tnp_0,lv4tnp_0,E_0,j_0,isss_0:=lmd_to_lv22tnp(lmd_0);
assert(isss_0);
E_0_4:=BaseChange(E_0,GF(p^4));

//take one basis of E_0[N_A]=(Z/N_A Z)^2.
P_A,Q_A:=ell_to_torsion_basis_2(E_0_4,N_A);
//take one basis of E_0[N_B]=(Z/N_B Z)^2.
P_B,Q_B:=ell_to_torsion_basis_2(E_0_4,N_B);
//========================================


//Bob calculates secretly.==================

coff_B:=Random({c: c in {0..(N_B-1)}});
R_B:=E_0_4!(P_B+coff_B*Q_B);
assert (Order(R_B) eq N_B);
E_B,PA_EB,QA_EB:=elliptic_isogeny_1ptker(E_0_4,R_B,P_A,Q_A);
assert(Order(PA_EB) eq N_A);
assert(Order(QA_EB) eq N_A);

//Note that the following data are public.
//E_A,PB_EA,QB_EA,E_B,PA_EB,QA_EB;
//===========================================




//construction auxiliary poinsts.================
"attack start.";

add_ct:=0;

a:=N_A-N_B; 
assert(a gt 0);



//if a is not squre.
E_pr,alpha_P_A,alpha_Q_A:=construct_auxiliary_img_6(E_0_4,N_A,N_B,P_A,Q_A);


/*
//if a is squre.
assert(IsSquare(a)); 
b:=IntegerRing()!Sqrt(a);  
alpha_0:=MultiplicationByMMap(E_0_4,b);  //4倍算 alpha_0:E_0->E_0.
E_pr:=E_0_4;
alpha_P_A:=alpha_0(P_A);
alpha_Q_A:=alpha_0(Q_A);
*/
//==============================================


atk_gen:=main_torsion_attack_3(E_0_4,E_B,E_pr,N_A,N_B,P_A,Q_A,PA_EB,QA_EB,alpha_P_A,alpha_Q_A,precomp_for_N_A);


Atk_Ker_phB:={k*E_0_4!atk_gen:k in {0..N_B}};  //Attacker.
Bob_Ker_phB:={k*E_0_4!R_B:k in {0..N_B}};      //Bob.
"";
"attack result.",Atk_Ker_phB eq Bob_Ker_phB;


//--------------------

