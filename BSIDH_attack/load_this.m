//please load this file.

load "setting.m";  
load "func_additions.m";
load "func_isogeny.m";
load "func_elliptic_theta.m";
load "func_theta_trans.m";
load "func_torsion_attack.m";

//--------------------------

//input p,N_A,N_B here.

p:=202546499;//27bit
N_A:=7^2 * 11^2 * 19 * 29 * 31;
N_B:=2^2 * 3 * 5^3 * 13^2 * 17;


load "test_BSIDH_attack.m"; 
