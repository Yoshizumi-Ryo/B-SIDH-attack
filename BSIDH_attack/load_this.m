//please load this file.

load "setting.m";  
load "func_additions.m";
load "func_isogeny.m";
load "func_elliptic_theta.m";
load "func_theta_trans.m";
load "func_torsion_attack.m";
//--------------------------

//input p,N_A,N_B here.
p:=911;
N_A:=5*7*13;
N_B:=2^4*3;

load "test_BSIDH_attack.m"; 
