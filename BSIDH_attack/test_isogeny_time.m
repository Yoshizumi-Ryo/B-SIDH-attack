//calculate time of l-isogeny of dim2.

load "setting.m";  

load "func_Mum_to_theta.m";
load "func_additions.m";
load "func_isogeny.m";
load "func_elliptic_theta.m";
load "func_theta_trans.m";
load "func_torsion_attack.m";


//calculate time of l-isogeny of the point with order N. 
procedure compute_isogeny(p,l)
  assert(IsPrime(l));
  assert(p mod 4 eq 3);
  assert((IsDivisibleBy(p+1,l)) or (IsDivisibleBy(p-1,l)));
  _<x>:=PolynomialRing(GF(p));
  E_0:=EllipticCurve(x^3-x);

  N:=1;
  for NN in {1..10} do
    if IsPrime(NN) then
      if NN ne l then
        N:=N*NN;
      end if;
    end if;
  end for;
  
  _<x>:=PolynomialRing(GF(p^4));
  E_0_4:=EllipticCurve(x^3-x);


  //"CI1";
  //N;
  P,Q:=ell_to_torsion_basis_2(E_0_4,N);
  //"CI2";


  lmd_0,lv22tnp_0,lv4tnp_0,E_0_4,j_0,isss_0:=E_to_lmd(E_0_4);
  S1,S2:=ell_to_torsion_basis_2(E_0_4,l);
  S12:=S1+S2;

  lv4tc_S1:=uvw_to_lv4tc(lmd_0,lv22tnp_0,S1[1],S1[2],S1[3]);
  lv4tc_S2:=uvw_to_lv4tc(lmd_0,lv22tnp_0,S2[1],S2[2],S2[3]);
  lv4tc_S12:=uvw_to_lv4tc(lmd_0,lv22tnp_0,S12[1],S12[2],S12[3]);
  lv4tc_P:=uvw_to_lv4tc(lmd_0,lv22tnp_0,P[1],P[2],P[3]);
  lv4tc_Q:=uvw_to_lv4tc(lmd_0,lv22tnp_0,Q[1],Q[2],Q[3]);
  lv4tc_PpS1:=uvw_to_lv4tc(lmd_0,lv22tnp_0,(P+S1)[1],(P+S1)[2],(P+S1)[3]);
  lv4tc_QpS1:=uvw_to_lv4tc(lmd_0,lv22tnp_0,(Q+S1)[1],(Q+S1)[2],(Q+S1)[3]);
  lv4tc_PpS2:=uvw_to_lv4tc(lmd_0,lv22tnp_0,(P+S2)[1],(P+S2)[2],(P+S2)[3]);
  lv4tc_QpS2:=uvw_to_lv4tc(lmd_0,lv22tnp_0,(Q+S2)[1],(Q+S2)[2],(Q+S2)[3]);


  lv4tc_e1:=ell_prod_lv4tc(lv4tc_S1,lv4tc_S1); //e1=(S1,S1)
  lv4tc_e2:=ell_prod_lv4tc(lv4tc_S2,lv4tc_S2); //e2=(S2,S2)
  lv4tc_e12:=ell_prod_lv4tc(lv4tc_S12,lv4tc_S12);
  lv4tnp_dm:=ell_prod_lv4tc(lv4tnp_0,lv4tnp_0);
  lv4tc_x:=ell_prod_lv4tc(lv4tc_P,lv4tc_Q);  //x=(P,Q).
  lv4tc_xpe1:=ell_prod_lv4tc(lv4tc_PpS1,lv4tc_QpS1); //x+e1=(P+S1,Q+S1).
  lv4tc_xpe2:=ell_prod_lv4tc(lv4tc_PpS2,lv4tc_QpS2); //x+e2=(P+S2,Q+S2).

  //"CI3";
  //---------------------------
  time_nullpt:=Time();
  time_nullpt_1:=Time();
  Mat_F:=const_Mat_F(l);
  r,set_vec_t,index_j:=const_index_t_j_3(l,Mat_F);
  "time_precomp.",Time(time_nullpt_1);

  time_nullpt_3:=Time();
  lv4tnp_cd:=lv4tnp_of_codomain(l,r,set_vec_t,index_j,lv4tnp_dm,lv4tc_e1,lv4tc_e2,lv4tc_e12);
  "time_isogeny.";Time(time_nullpt_3);
  "1.time_null_point.",Time(time_nullpt);
  //---------------------------
  
  //-----------------------------
  time_non_nullpt_A:=Time();
  time_non_nullpt_B:=Time();
  lincom_e1e2:=linear_combination(lv4tnp_dm,l,lv4tc_e1,lv4tc_e2,lv4tc_e12); 
  "2.time_for_all_pts.",Time(time_non_nullpt_B);
  
  time_non_nullpt_C:=Time();
  lv4tnp_imgx:=image_of_point(lincom_e1e2,l,Mat_F,set_vec_t,index_j,lv4tnp_dm,lv4tc_e1,lv4tc_e2,lv4tc_e12,lv4tc_x,lv4tc_xpe1,lv4tc_xpe2);
  "3.time_for_the_point",Time(time_non_nullpt_C);
  "4.time_non_null_point.",Time(time_non_nullpt_A);
  //-----------------------------
end procedure;



//===================
//example. 

p:=18628989148679788872005065350440589045599;

p:=826791736418446924644415105270960270928927659729776400179861442336062222833458285859;

//----------------------

for l in {3..50} do
  if IsPrime(l) then
    assert(IsDivisibleBy(p+1,l) or IsDivisibleBy(p-1,l));
    l;
    compute_isogeny(p,l);
    "";
  end if;
end for;








//---------------

/*
実装結果の見方について.

1.time_null_point.     :codomainのtheta null pointを計算するのにかかった時間(秒).

4.time_non_null_point. :null pointでない元xの像を計算するのにかかった時間(秒).

上の時間4は次の2つの時間の合計になっています.

2.time_for_all_pts.    :xに依存しない計算にかかった時間(秒).i.e. 0でないx_1,x_2,...に対して, ここの計算は一度だけすれば十分です.

3.time_for_the_point   :xに依存する計算にかかった時間(秒)

先日送ったpdfのtableにおいて, null ptは上の時間1, not-null ptは時間3を意味しています.
*/

