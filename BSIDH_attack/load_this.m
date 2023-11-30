//please load this file.

//input p characteristic.

/*
p:=202546499;
N_A:=7^2 * 11^2 * 19 * 29 * 31;
N_B:=2^2 * 3 * 5^3 * 13^2 * 17;
*/


p:=202546499;
N_A:=7^2 * 11^2 * 19 * 29 * 31;
N_B:=2^2 * 3 * 5^3 * 13^2 * 17;


load "setting.m";  
load "func_Mum_to_theta.m";
load "func_additions.m";
load "func_isogeny.m";
load "func_elliptic_theta.m";
load "func_theta_trans.m";
load "func_torsion_attack.m";

//input N_A,N_B here.

load "test_BSIDH_attack.m"; // other parameters here.
