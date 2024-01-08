//test space for [LR22].

load "setting.m";  

load "func_Mum_to_theta.m";
load "func_additions.m";
load "func_isogeny.m";
load "func_elliptic_theta.m";
load "func_theta_trans.m";
load "func_torsion_attack.m";
load "func_LRtwotwo.m";



//20ビット安全.



p:=10259;
N_A:=19;
N_B:=15;
l:=19;

l mod 4;
assert(IsDivisibleBy(N_A,l));

load "const_isotropic.m"; 

//lv4tnp,lv4_e1,lv4_e2,lv4_e12,lv4_x,lv4_xpe1,lv4_xpe2,l:=



procedure base_field(lv4tc)
  for key in lv4keys do
    Parent(lv4_xpe1[key]);
  end for;
end procedure;




assert(Is_lv4tnp(lv4tnp));
assert(IsOrder(lv4tnp,lv4_e1,l));
assert(IsOrder(lv4tnp,lv4_e2,l));
assert(IsOrder(lv4tnp,lv4_e12,l));
assert(IsOrder(lv4tnp,lv4_x,N_B));


//null point.
//-----------
//way of [CR15].
time0:=Time();
Mat_F:=const_Mat_F(l);
r,set_vec_t,index_j:=const_index_t_j_3(l,Mat_F);
lv4tnp_cod_1:=lv4tnp_of_codomain(l,r,set_vec_t,index_j,lv4tnp,lv4_e1,lv4_e2,lv4_e12);
assert(Is_lv4tnp(lv4tnp_cod_1));
Time(time0);

//------------
//way of [LR22]
time1:=Time();
lv4tnp_cod_2,lv4_e1,lv4_e2,lv4_e12:=codomain_tnp(lv4tnp,lv4_e1,lv4_e2,lv4_e12,l);
Time(time1);
assert(Is_lv4tnp(lv4tnp_cod_2));
lv22tnp_cod_2:=to_lv22(lv4tnp_cod_2);

//----------
//way of [LR22]new

time1:=Time();
lv4tnp_cod_3,mu_1_lpow,mu_2_lpow,mu_12_lpow:=codomain_tnp_2(lv4tnp,lv4_e1,lv4_e2,lv4_e12,l);
Time(time1);
assert(Is_lv4tnp(lv4tnp_cod_3));
//-----------

//-------------
assert(eq_tc(lv4tnp_cod_1,lv4tnp_cod_2));
assert(eq_tc(lv4tnp_cod_2,lv4tnp_cod_3));
assert(eq_tc(lv4tnp_cod_1,lv4tnp_cod_3));
//-------------





//for general point.-------------
//[CR15].
time0:=Time();
lincom_e1e2:=linear_combination(lv4tnp,l,lv4_e1,lv4_e2,lv4_e12); 
lv4_fx_1:=image_of_point(lincom_e1e2,l,Mat_F,set_vec_t,index_j,lv4tnp,lv4_e1,lv4_e2,lv4_e12,lv4_x,lv4_xpe1,lv4_xpe2);
Time(time0);
assert(IsOrder(lv4tnp_cod_2,lv4_fx_1,N_B));
//--------------
//[LR22].
time0:=Time();
lv4_fx_2:=ll_isogeny(lv4tnp,lv4_e1,lv4_e2,lv4_e12,lv4_x,lv4_xpe1,lv4_xpe2,l);
Time(time0);
assert(IsOrder(lv4tnp_cod_2,lv4_fx_2,N_B));
//--------------
//[LR22]new.
time0:=Time();
lv4_fx_3:=ll_isogeny_2(lv4tnp,lv4_e1,lv4_e2,lv4_e12,lv4_x,lv4_xpe1,lv4_xpe2,l,mu_1_lpow,mu_2_lpow,mu_12_lpow);
Time(time0);


//----------
eq_tc(lv4_fx_1,lv4_fx_2) or eq_tc(lv4_fx_1,const_inv_lv4tc(lv4_fx_2));
eq_tc(lv4_fx_2,lv4_fx_3);
//--------------






//-------------------
//(Z/4Z)^2の自己同型について. 

matset:={};
for abcd in CartesianPower({0..3},4) do
  a:=abcd[1];
  b:=abcd[2];
  c:=abcd[3];
  d:=abcd[4];
  if ((a*d-b*c) mod 2) eq 1 then
    abcd;
    matset join:={abcd};
  end if;
end for;




for abcd in matset do
  a:=abcd[1];
  b:=abcd[2];
  c:=abcd[3];
  d:=abcd[4];
  new_lv4_fx_2:=AssociativeArray();
  for key in lv4keys do
    matkey:=[(a*key[1]+b*key[2]) mod 4,(c*key[1]+d*key[2]) mod 4];
    new_lv4_fx_2[matkey]:=lv4_fx_2[key];
  end for;
  if eq_tc(lv4_fx_1,new_lv4_fx_2) then
    a,b,c,d;
  end if;
  assert (Keys(new_lv4_f0) eq lv4keys);
end for;
//----------------




_<x>:=PolynomialRing(GF(p));
E_001:=EllipticCurve(x^3-x);
E_00t1:=EllipticCurve(x^3+x);
#E_001 eq (p+1);
#E_00t1 eq (p+1);
IsIsomorphic(E_001,E_00t1);


E_002:=BaseChange(E_001,GF(p^2));
E_00t2:=BaseChange(E_00t1,GF(p^2));
#E_002 eq (p+1)^2;
#E_00t2 eq (p+1)^2;
IsIsomorphic(E_002,E_00t2);


_<t>:=PolynomialRing(GF(p^2));
E_01:=EllipticCurve(t^3+t);
QuadraticTwists(E_01);



