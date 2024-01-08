//please load this file for B-SIDH attack.

load "setting.m";  

load "func_additions.m";
load "func_isogeny.m";
load "func_elliptic_theta.m";
load "func_theta_trans.m";
load "func_torsion_attack.m";

//--------------------------

//input p,N_A,N_B here.

//30bit security.
p := 276154505650672190920223;
N_A:=13 * 23 * 37 * 43 * 47^2 * 59 * 61 * 71 * 73 * 97 * 101;
N_B:=2^5 * 3^6 * 7 * 11^4 * 19 * 29^3 * 67 * 79;

load "test_BSIDH_attack.m"; 
