//=======================
//======================

//test page of Kummer line and Kummer surface.

load "setting.m";  
load "func_Mum_to_theta.m";
load "func_additions.m";
load "func_isogeny.m";
load "func_elliptic_theta.m";
load "func_theta_trans.m";
load "func_torsion_attack.m";


p:=625750366823999;
//here,input p.

load "func_Kummer.m";


iso_lmd_set:=supersinular_Legendre(p);

lm1:=Random(Random(iso_lmd_set));
lm2:=Random(Random(iso_lmd_set));


_<x>:=PolynomialRing(GF(p^2));

lv2tnp1,sqrt_lm1,sqrt_lmm11:=LegendreEll_to_lv2tnp(lm1);
lv2tnp2,sqrt_lm2,sqrt_lmm12:=LegendreEll_to_lv2tnp(lm2);

E1:=EllipticCurve(x*(x-1)*(x-lm1));
E2:=EllipticCurve(x*(x-1)*(x-lm2));

assert(IsSupersingular(E1));
assert(IsSupersingular(E2));


l:=11;
assert(IsDivisibleBy(p+1,l));
E1_e1,E1_e2:=ell_to_torsion_basis_2(E1,l);
E2_e1,E2_e2:=ell_to_torsion_basis_2(E2,l);


E1_x:=Random(E1);
E2_x:=Random(E2);

//-------------------
//dim1 theta.

lv2tc_E1_e1   :=Legendre_to_lv2tc(lv2tnp1,sqrt_lm1,sqrt_lmm11,E1_e1);
lv2tc_E1_e2   :=Legendre_to_lv2tc(lv2tnp1,sqrt_lm1,sqrt_lmm11,E1_e2);
lv2tc_E1_e12  :=Legendre_to_lv2tc(lv2tnp1,sqrt_lm1,sqrt_lmm11,(E1_e1+E1_e2));
lv2tc_E1_x    :=Legendre_to_lv2tc(lv2tnp1,sqrt_lm1,sqrt_lmm11,(E1_x));
lv2tc_E1_xpe1 :=Legendre_to_lv2tc(lv2tnp1,sqrt_lm1,sqrt_lmm11,(E1_x+E1_e1));
lv2tc_E1_xpe2 :=Legendre_to_lv2tc(lv2tnp1,sqrt_lm1,sqrt_lmm11,(E1_x+E1_e2));
lv2tc_E1_2x   :=Legendre_to_lv2tc(lv2tnp1,sqrt_lm1,sqrt_lmm11,(2*E1_x));
lv2tc_E1_2xpe1:=Legendre_to_lv2tc(lv2tnp1,sqrt_lm1,sqrt_lmm11,(2*E1_x+E1_e1));
lv2tc_E1_2xpe2:=Legendre_to_lv2tc(lv2tnp1,sqrt_lm1,sqrt_lmm11,(2*E1_x+E1_e2));

lv2tc_E2_e1   :=Legendre_to_lv2tc(lv2tnp2,sqrt_lm2,sqrt_lmm12,E2_e1);
lv2tc_E2_e2   :=Legendre_to_lv2tc(lv2tnp2,sqrt_lm2,sqrt_lmm12,E2_e2);
lv2tc_E2_e12  :=Legendre_to_lv2tc(lv2tnp2,sqrt_lm2,sqrt_lmm12,(E2_e1+E2_e2));
lv2tc_E2_x    :=Legendre_to_lv2tc(lv2tnp2,sqrt_lm2,sqrt_lmm12,(E2_x));
lv2tc_E2_xpe1 :=Legendre_to_lv2tc(lv2tnp2,sqrt_lm2,sqrt_lmm12,(E2_x+E2_e1));
lv2tc_E2_xpe2 :=Legendre_to_lv2tc(lv2tnp2,sqrt_lm2,sqrt_lmm12,(E2_x+E2_e2));
lv2tc_E2_2x   :=Legendre_to_lv2tc(lv2tnp2,sqrt_lm2,sqrt_lmm12,(2*E2_x));
lv2tc_E2_2xpe1:=Legendre_to_lv2tc(lv2tnp2,sqrt_lm2,sqrt_lmm12,(2*E2_x+E2_e1));
lv2tc_E2_2xpe2:=Legendre_to_lv2tc(lv2tnp2,sqrt_lm2,sqrt_lmm12,(2*E2_x+E2_e2));


//----------------
//2dim.



lv2tnp  :=product_theta(lv2tnp1      ,lv2tnp2);
lv2_e1  :=product_theta(lv2tc_E1_e1  ,lv2tnp2);
lv2_e2  :=product_theta(lv2tnp1      ,lv2tc_E2_e2);
lv2_e12 :=product_theta(lv2tc_E1_e1  ,lv2tc_E2_e2);
lv2_x   :=product_theta(lv2tc_E1_x   ,lv2tc_E2_x);
lv2_xpe1:=product_theta(lv2tc_E1_xpe1,lv2tc_E2_x);
lv2_xpe2:=product_theta(lv2tc_E1_x   ,lv2tc_E2_xpe2);

lv2_2x   :=product_theta(lv2tc_E1_2x   ,lv2tc_E2_2x);
lv2_2xpe1:=product_theta(lv2tc_E1_2xpe1,lv2tc_E2_2x);
lv2_2xpe2:=product_theta(lv2tc_E1_2x   ,lv2tc_E2_2xpe2);



//some checking.

lv2_xpe12:=Extended_Add(lv2tnp,lv2_x,lv2_e1,lv2_e2,lv2_xpe1,lv2_e12,lv2_xpe2);

assert(eq_lv2tc(lv2tnp,mult_g2_n2(lv2tnp,l,lv2_e1)));
assert(eq_lv2tc(lv2tnp,mult_g2_n2(lv2tnp,l,lv2_e1)));
assert(eq_lv2tc(lv2tnp,mult_g2_n2(lv2tnp,l,lv2_e12)));

assert(eq_lv2tc(lv2_x,xpke(lv2tnp,l,lv2_x,lv2_e1,lv2_xpe1)));
assert(eq_lv2tc(lv2_x,xpke(lv2tnp,l,lv2_x,lv2_e2,lv2_xpe2)));
assert(eq_lv2tc(lv2_x,xpke(lv2tnp,l,lv2_x,lv2_e12,lv2_xpe12)));

lv2_2x_2:=mult_g2_n2(lv2tnp,2,lv2_x);
assert(eq_lv2tc(lv2_2x,lv2_2x_2));


assert(eq_Assoc(Ladder_g2n2(lv2tnp,7,lv2_x),mult_g2_n2(lv2tnp,7,lv2_x)));
assert(eq_Assoc(xpke(lv2tnp,15,lv2_x,lv2_e1,lv2_xpe1),xpke_2(lv2tnp,lv2_x,15,lv2_e1,lv2_xpe1)));





//---------------------
//isogeny of tnp.

lv2tnp_cd,mu1_lpow,mu2_lpow,mu12_lpow:=codomain_lv2tnp(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l);

//check.
lv2tnp_cd_2:=DiffAdd_lv2_linsys_2(lv2tnp_cd,lv2tnp_cd,lv2tnp_cd,lv2tnp_cd);
eq_lv2tc(lv2tnp_cd,lv2tnp_cd_2);

//------------------
//image of general point.

lv2_f0:=image_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,lv2tnp,lv2_e1,lv2_e2,mu1_lpow,mu2_lpow,mu12_lpow);
assert(eq_lv2tc(lv2tnp_cd,lv2_f0));

lv2_fx:=image_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,lv2_x,lv2_xpe1,lv2_xpe2,mu1_lpow,mu2_lpow,mu12_lpow);
get_order_lv2(lv2tnp,lv2_x);
get_order_lv2(lv2tnp_cd,lv2_fx);

lv2_f2x:=image_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,lv2_2x,lv2_2xpe1,lv2_2xpe2,mu1_lpow,mu2_lpow,mu12_lpow);
lv2_2fx:=mult_g2_n2(lv2tnp_cd,2,lv2_fx);
assert(eq_lv2tc(lv2_f2x,lv2_2fx));

//=========================

