
load "setting.m";  
load "func_additions.m";
load "func_isogeny.m";
load "func_elliptic_theta.m";
load "func_theta_trans.m";
load "func_torsion_attack.m";

//------------------------


//30bit security.
p := 276154505650672190920223;
N_A:=13 * 23 * 37 * 43 * 47^2 * 59 * 61 * 71 * 73 * 97 * 101;
N_B:=2^5 * 3^6 * 7 * 11^4 * 19 * 29^3 * 67 * 79;


load "func_Kummer.m";
load "func_attack_Kum.m";
load "test_BSIDH_attack_Kum.m"; 

