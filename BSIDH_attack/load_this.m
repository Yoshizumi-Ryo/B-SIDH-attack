//please load this file.

//input p characteristic.


p:=1479871;
N_A:=3^6*5*7*29;
N_B:=1216;

load "setting.m";  
load "func_Mum_to_theta.m";
load "func_additions.m";
load "func_isogeny.m";
load "func_elliptic_theta.m";
load "func_theta_trans.m";
load "func_torsion_attack.m";

//input N_A,N_B here.

load "test_BSIDH_attack.m"; // other parameters here.
