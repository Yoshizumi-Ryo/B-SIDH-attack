//==============================================
//start of torsion_attack6.m


//auxi. functions.----------------------------


function Hadamard_transform(lv2tc)
  t_lv2tc:=AssociativeArray();
  x:=lv2tc[[0,0]];
  y:=lv2tc[[1,0]];
  z:=lv2tc[[0,1]];
  w:=lv2tc[[1,1]];
  t_lv2tc[[0,0]]:=x+y+z+w;
  t_lv2tc[[1,0]]:=x-y+z-w;
  t_lv2tc[[0,1]]:=x+y-z-w;
  t_lv2tc[[1,1]]:=x-y-z+w;
  return t_lv2tc;
end function;




//最後に楕円直積をthetaの積に変換する.
//cf.[DMPR23].p18.
function theta_product_lv2(lv2tnp,lv2_x,lv2_y,zeta_4)
  assert(Is_prod_ell(lv2tnp_to_lv22tnpsq(lv2tnp)));
  //check.
  if Is_tehta_prod_ell(lv2tnp_to_lv22tnpsq(lv2tnp)) then
    return lv2tnp,lv2_x,lv2_y;
  else
    _,evenkey:=Is_tehta_prod_ell(lv2tnp_to_lv22tnpsq(lv2tnp));
  end if;
  //setp1.
  //"before step1.",evenkey;
  if evenkey eq [0,0,0,0] then
    lv2tnp[[1,0]]*:=zeta_4;
    lv2tnp[[0,1]]*:=zeta_4;
    lv2_x[[1,0]]*:=zeta_4;
    lv2_x[[0,1]]*:=zeta_4;
    lv2_y[[1,0]]*:=zeta_4;
    lv2_y[[0,1]]*:=zeta_4;
  end if;
  //cehck.
  if Is_tehta_prod_ell(lv2tnp_to_lv22tnpsq(lv2tnp)) then
    return lv2tnp,lv2_x,lv2_y;
  else
    _,evenkey:=Is_tehta_prod_ell(lv2tnp_to_lv22tnpsq(lv2tnp));
  end if;
  //"before step2.",evenkey;
  //step2.
  if (evenkey[3] eq 0) and  (evenkey[4] eq 0) then
    assert((evenkey[1] ne 0) or (evenkey[2] ne 0));
    lv2tnp:=Hadamard_transform(lv2tnp);
    lv2_x:=Hadamard_transform(lv2_x);
    lv2_y:=Hadamard_transform(lv2_y);
  end if;
  //cehck.
  if Is_tehta_prod_ell(lv2tnp_to_lv22tnpsq(lv2tnp)) then
    return lv2tnp,lv2_x,lv2_y;
  else
    _,evenkey:=Is_tehta_prod_ell(lv2tnp_to_lv22tnpsq(lv2tnp));
  end if;
  assert((evenkey[3] ne 0) or  (evenkey[4] ne 0));
  //step3.
  //"before step3.",evenkey;
  if ((evenkey[3] eq 0) and  (evenkey[4] eq 1)) then
    resotor:=lv2tnp[[0,1]];
    lv2tnp[[0,1]]:=lv2tnp[[1,1]];
    lv2tnp[[1,1]]:=resotor;
    resotor:=lv2_x[[0,1]];
    lv2_x[[0,1]]:=lv2_x[[1,1]];
    lv2_x[[1,1]]:=resotor;
    resotor:=lv2_y[[0,1]];
    lv2_y[[0,1]]:=lv2_y[[1,1]];
    lv2_y[[1,1]]:=resotor;
  elif ((evenkey[3] eq 1) and  (evenkey[4] eq 0)) then
    resotor:=lv2tnp[[1,0]];
    lv2tnp[[1,0]]:=lv2tnp[[1,1]];
    lv2tnp[[1,1]]:=resotor;
    resotor:=lv2_x[[1,0]];
    lv2_x[[1,0]]:=lv2_x[[1,1]];
    lv2_x[[1,1]]:=resotor;
    resotor:=lv2_y[[1,0]];
    lv2_y[[1,0]]:=lv2_y[[1,1]];
    lv2_y[[1,1]]:=resotor;
  end if;
  //cehck.
  if Is_tehta_prod_ell(lv2tnp_to_lv22tnpsq(lv2tnp)) then
    return lv2tnp,lv2_x,lv2_y;
  else
    _,evenkey:=Is_tehta_prod_ell(lv2tnp_to_lv22tnpsq(lv2tnp));
  end if;
  //"before step4.",evenkey;
  assert(evenkey eq [0,0,1,1]);
  //step4.
  lv2tnp[[1,0]]*:=zeta_4;
  lv2tnp[[0,1]]*:=zeta_4;
  lv2_x[[1,0]]*:=zeta_4;
  lv2_x[[0,1]]*:=zeta_4;
  lv2_y[[1,0]]*:=zeta_4;
  lv2_y[[0,1]]*:=zeta_4;
  assert(Is_tehta_prod_ell(lv2tnp_to_lv22tnpsq(lv2tnp)));
  return lv2tnp,lv2_x,lv2_y;
end function;






function prod_lv2tc(E1,E2,S1,S2,lv2tnp_1,sqrt_lm_1,sqrt_lmm1_1,lv2tnp_2,sqrt_lm_2,sqrt_lmm1_2)
  assert(S1 in E1);
  assert(S2 in E2);
  //assert(S1[1] in GF(p^2));
  //assert(S2[1] in GF(p^2));
  lv2_S1:=Legendre_to_lv2tc(lv2tnp_1,sqrt_lm_1,sqrt_lmm1_1,S1);
  lv2_S2:=Legendre_to_lv2tc(lv2tnp_2,sqrt_lm_2,sqrt_lmm1_2,S2);
  lv2_S:=product_theta(lv2_S1,lv2_S2);
  return lv2_S;
end function;









//aux alpha:E_0_4->E_pr.
function MainTorsionAttackKummer(E_0_4,E_B,E_pr,N_A,N_B,P_A,Q_A,PA_EB,QA_EB,alpha_P_A,alpha_Q_A,zeta_8)
  "p=      ",p;
  "N_A=    ",N_A;
  "N_B=    ",N_B;
  zeta_4:=zeta_8^2;
  time_total:=Time();
  time_prepare:=Time();
  assert(N_A gt N_B);
  assert(P_A in E_0_4);
  assert(Q_A in E_0_4);
  assert(PA_EB in E_B);
  assert(QA_EB in E_B);
  assert(alpha_P_A in E_pr);
  assert(alpha_Q_A in E_pr);
  lmd_0,lv22tnp_0,lv4tnp_0,E_0_4,j_0,isss_0:=E_to_lmd(E_0_4);
  //"wait for chaing elliptic curve to another isomorphic one."; 
  time_ell_change:=Time();
  lmd_B,lv22tnp_B,lv4tnp_B,E_B2,j_B,isss_B,iso_E_B_2:=E_to_lmd(E_B);
  lmd_pr,lv22tnp_pr,lv4tnp_pr,E_pr,j_pr,isss_pr,iso_E_pr:=E_to_lmd(E_pr);
  //"fin.",Time(time_ell_change);
  assert(isss_0);
  assert(isss_B);
  assert(isss_pr);

  //(N_A,N_A)-isogeny F:E_cd*E_B->E_0*E'_B.
  //basis_KerF:={[alpha(P_A),PA_EB],[alpha(Q_A),QA_EB]}; //in E_cd*E_B.
  //we will call e_1=[alpha(P_A),PA_EB], e_2=[alpha(Q_A),QA_EB].
  //Next we want to calculate  F(0,PA_EB)=(S_1,*), F(0,QA_EB)=(S_2,*), because Ker(phi_B)=<S_1,S_2>.
  S1,S2:=ell_to_torsion_basis_2(E_B2,N_B); //attacker will use.

  _<x>:=PolynomialRing(GF(p^2));
  //E_B_2ex:=EllipticCurve(x*(x-1)*(x-GF(p^2)!lmd_B));
  
  //==========================
  //Theta.
  lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr:=LegendreEll_to_lv2tnp(lmd_pr);
  lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B :=LegendreEll_to_lv2tnp(lmd_B);
  lv2tnp_EprEB:=product_theta(lv2tnp_Epr,lv2tnp_EB);
  lv2tnp_EprEB_2:=prod_lv2tc(E_pr,E_B2,E_pr!0,E_B2!0,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);
  
  f1_Epr:=iso_E_pr(alpha_P_A);    
  f1_EB:=iso_E_B_2(PA_EB); 
  lv2_f1:=prod_lv2tc(E_pr,E_B2,f1_Epr,f1_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);
  
  f2_Epr:=iso_E_pr(alpha_Q_A);    
  f2_EB:=iso_E_B_2(QA_EB);    
  lv2_f2:=prod_lv2tc(E_pr,E_B2,f2_Epr,f2_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);

  f12_Epr:=f1_Epr+f2_Epr;
  f12_EB:=f1_EB+f2_EB;
  lv2_f12:=prod_lv2tc(E_pr,E_B2,f12_Epr,f12_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);

  //------------ZS1.
  ZS1_Epr:=E_pr!0;
  ZS1_EB:=S1;
  lv2_ZS1:=prod_lv2tc(E_pr,E_B2,ZS1_Epr,ZS1_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);

  ZS1f1_Epr:=f1_Epr;
  ZS1f1_EB:=f1_EB+S1;
  lv2_ZS1pf1:=prod_lv2tc(E_pr,E_B2,ZS1f1_Epr,ZS1f1_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);

  ZS1f2_Epr:=f2_Epr;
  ZS1f2_EB:=f2_EB+S1;
  lv2_ZS1pf2:=prod_lv2tc(E_pr,E_B2,ZS1f2_Epr,ZS1f2_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);
  //--------------

  //------------ZS2.
  ZS2_Epr:=E_pr!0;
  ZS2_EB:=S2;
  lv2_ZS2:=prod_lv2tc(E_pr,E_B2,ZS2_Epr,ZS2_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);

  ZS2f1_Epr:=f1_Epr;
  ZS2f1_EB:=f1_EB+S2;
  lv2_ZS2pf1:=prod_lv2tc(E_pr,E_B2,ZS2f1_Epr,ZS2f1_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);

  ZS2f2_Epr:=f2_Epr;
  ZS2f2_EB:=f2_EB+S2;
  lv2_ZS2pf2:=prod_lv2tc(E_pr,E_B2,ZS2f2_Epr,ZS2f2_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);
  //===========================


  //---------------------------
  lv2_img0      :=lv2tnp_EprEB;
  lv2_img_f1    :=lv2_f1;
  lv2_img_f2    :=lv2_f2;
  lv2_img_f12   :=lv2_f12;
  lv2_img_ZS1   :=lv2_ZS1;
  lv2_img_ZS1pf1:=lv2_ZS1pf1;
  lv2_img_ZS1pf2:=lv2_ZS1pf2;
  lv2_img_ZS2   :=lv2_ZS2;
  lv2_img_ZS2pf1:=lv2_ZS2pf1;
  lv2_img_ZS2pf2:=lv2_ZS2pf2;
  //-----------------------------
 

  "E_0:lmd=",lmd_0;
  "E' :lmd=",lmd_pr;//'
  "E_B:lmd=",lmd_B;
  "";
  "Time for prepare",Time(time_prepare);
  //=========================================
  time_isogeny:=Time();
  "start attack isogeny-------------------------";
  fac:=decomposition_to_seq(N_A); //the order of composition.
  fac;
  s:=1;

  for i in {1..#fac} do
    //"";
    time_domain:=Time();
    l:=fac[i];
    kk:=IntegerRing()!(N_A/(s*l));

    //"we calculated degree.", cut_out_first(fac,(i-1));
    //"from now, we calculate degree l=",l;
    "l=",l;
    //"we remain degree.", cut_out_last(fac,(#fac-i));

    assert(s*l*kk eq N_A);
  
    //-----------------------------
    lv2tnp    :=lv2_img0;
    lv2_f1    :=lv2_img_f1;
    lv2_f2    :=lv2_img_f2;
    lv2_f12   :=lv2_img_f12;
    lv2_ZS1   :=lv2_img_ZS1;
    lv2_ZS1pf1:=lv2_img_ZS1pf1;
    lv2_ZS1pf2:=lv2_img_ZS1pf2;
    lv2_ZS2   :=lv2_img_ZS2;
    lv2_ZS2pf1:=lv2_img_ZS2pf1;
    lv2_ZS2pf2:=lv2_img_ZS2pf2;
    //-----------------------------

    
    //"order_check1.",IsOrder_g2lv2(lv4tnp_dm,N_B,lv4tc_0S1_dm);
    //"order_check2.",IsOrder_g2lv2(lv4tnp_dm,l*kk,lv4tc_f1_dm);


    //"wait for calculate on the domain.";
    time_A:=Time();
    //"MT15";
    lincom_f1f2:=AssociativeArray();
    lincom_f1f2[[1,1]]:=lv2_f12;
    //"MT15.1";
    lincom_f1f2[[kk,0]]    :=mult_g2_n2(lv2tnp,kk  ,lv2_f1);
    lincom_f1f2[[(kk+1),0]]:=mult_g2_n2(lv2tnp,kk+1,lv2_f1);
    lincom_f1f2[[0,kk]]    :=mult_g2_n2(lv2tnp,kk  ,lv2_f2);
    lincom_f1f2[[0,(kk+1)]]:=mult_g2_n2(lv2tnp,kk+1,lv2_f2);
    lincom_f1f2[[kk,kk]]   :=mult_g2_n2(lv2tnp,kk  ,lv2_f12);
    lincom_f1f2[[1,kk]]    :=xpke      (lv2tnp,kk  ,lv2_f1,lv2_f2,lv2_f12);
    lincom_f1f2[[1,kk+1]]  :=xpke      (lv2tnp,kk+1,lv2_f1,lv2_f2,lv2_f12);
    lincom_f1f2[[kk,1]]    :=xpke      (lv2tnp,kk  ,lv2_f2,lv2_f1,lv2_f12);
    lincom_f1f2[[kk+1,1]]  :=xpke      (lv2tnp,kk+1,lv2_f2,lv2_f1,lv2_f12);

    lv2_e1 :=lincom_f1f2[[kk,0]];
    lv2_e2 :=lincom_f1f2[[0,kk]];
    lv2_e12:=lincom_f1f2[[kk,kk]];

    //S1----------------------------------
    //"wait for calculate about S1.";
    time_S1:=Time(); 
    lv2_ZS1pf12:=Extended_Add(lv2tnp,lv2_ZS1,lv2_f1,lv2_f2,lv2_ZS1pf1,lv2_f12,lv2_ZS1pf2);
    tc_ZS1_lincomf1f2:=AssociativeArray();
    tc_ZS1_lincomf1f2[[0,0]]   :=lv2_ZS1;
    tc_ZS1_lincomf1f2[[1,0]]   :=lv2_ZS1pf1;
    tc_ZS1_lincomf1f2[[kk,0]]  :=xpke(lv2tnp,kk  ,lv2_ZS1   ,lv2_f1,lv2_ZS1pf1);
    tc_ZS1_lincomf1f2[[kk+1,0]]:=xpke(lv2tnp,kk+1,lv2_ZS1   ,lv2_f1,lv2_ZS1pf1);
    tc_ZS1_lincomf1f2[[0,1]]   :=lv2_ZS1pf2;
    tc_ZS1_lincomf1f2[[0,kk]]  :=xpke(lv2tnp,kk  ,lv2_ZS1   ,lv2_f2,lv2_ZS1pf2);
    tc_ZS1_lincomf1f2[[0,kk+1]]:=xpke(lv2tnp,kk+1,lv2_ZS1   ,lv2_f2,lv2_ZS1pf2);
    tc_ZS1_lincomf1f2[[kk,1]]  :=xpke(lv2tnp,kk  ,lv2_ZS1pf2,lv2_f1,lv2_ZS1pf12);
    tc_ZS1_lincomf1f2[[1,kk]]  :=xpke(lv2tnp,kk  ,lv2_ZS1pf1,lv2_f2,lv2_ZS1pf12);
    //S2----------------------------------
    //"wait for calculate about S2.";
    time_S2:=Time(); 
    lv2_ZS2pf12:=Extended_Add(lv2tnp,lv2_ZS2,lv2_f1,lv2_f2,lv2_ZS2pf1,lv2_f12,lv2_ZS2pf2);
    tc_ZS2_lincomf1f2:=AssociativeArray();
    tc_ZS2_lincomf1f2[[0,0]]   :=lv2_ZS2;
    tc_ZS2_lincomf1f2[[1,0]]   :=lv2_ZS2pf1;
    tc_ZS2_lincomf1f2[[kk,0]]  :=xpke(lv2tnp,kk  ,lv2_ZS2   ,lv2_f1,lv2_ZS2pf1);
    tc_ZS2_lincomf1f2[[kk+1,0]]:=xpke(lv2tnp,kk+1,lv2_ZS2   ,lv2_f1,lv2_ZS2pf1);
    tc_ZS2_lincomf1f2[[0,1]]   :=lv2_ZS2pf2;
    tc_ZS2_lincomf1f2[[0,kk]]  :=xpke(lv2tnp,kk  ,lv2_ZS2   ,lv2_f2,lv2_ZS2pf2);
    tc_ZS2_lincomf1f2[[0,kk+1]]:=xpke(lv2tnp,kk+1,lv2_ZS2   ,lv2_f2,lv2_ZS2pf2);
    tc_ZS2_lincomf1f2[[kk,1]]  :=xpke(lv2tnp,kk  ,lv2_ZS2pf2,lv2_f1,lv2_ZS2pf12);
    tc_ZS2_lincomf1f2[[1,kk]]  :=xpke(lv2tnp,kk  ,lv2_ZS2pf1,lv2_f2,lv2_ZS2pf12);
    //================================

    //img of 0.
    lv2_img0,mu1_lpow,mu2_lpow,mu12_lpow:=codomain_lv2tnp(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l);
    //---------
    //img of f_1.
    lv2_img_f1:=image_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,lv2_f1,lincom_f1f2[[(kk+1),0]],lincom_f1f2[[1,kk]],mu1_lpow,mu2_lpow,mu12_lpow);
    //img of f_2.
    lv2_img_f2:=image_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,lv2_f2,lincom_f1f2[[kk,1]],lincom_f1f2[[0,(kk+1)]],mu1_lpow,mu2_lpow,mu12_lpow);
    //img of f_12.
    lv2_img_f12:=image_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,lv2_f12,lincom_f1f2[[kk+1,1]],lincom_f1f2[[1,(kk+1)]],mu1_lpow,mu2_lpow,mu12_lpow);
    //--------------
    //img of ZS1.
    lv2_img_ZS1:=image_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,tc_ZS1_lincomf1f2[[0,0]],tc_ZS1_lincomf1f2[[kk,0]],tc_ZS1_lincomf1f2[[0,kk]],mu1_lpow,mu2_lpow,mu12_lpow);
    //img of ZS1+f_1.
    lv2_img_ZS1pf1:=image_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,tc_ZS1_lincomf1f2[[1,0]],tc_ZS1_lincomf1f2[[(kk+1),0]],tc_ZS1_lincomf1f2[[1,kk]],mu1_lpow,mu2_lpow,mu12_lpow);
    //img of ZS1+f_2.
    lv2_img_ZS1pf2:=image_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,tc_ZS1_lincomf1f2[[0,1]],tc_ZS1_lincomf1f2[[kk,1]],tc_ZS1_lincomf1f2[[0,(kk+1)]],mu1_lpow,mu2_lpow,mu12_lpow);
    //--------------
    //img of ZS2.
    lv2_img_ZS2:=image_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,tc_ZS2_lincomf1f2[[0,0]],tc_ZS2_lincomf1f2[[kk,0]],tc_ZS2_lincomf1f2[[0,kk]],mu1_lpow,mu2_lpow,mu12_lpow);
    //img of ZS2+f_1.
    lv2_img_ZS2pf1:=image_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,tc_ZS2_lincomf1f2[[1,0]],tc_ZS2_lincomf1f2[[(kk+1),0]],tc_ZS2_lincomf1f2[[1,kk]],mu1_lpow,mu2_lpow,mu12_lpow);
    //img of ZS2+f_2.
    lv2_img_ZS2pf2:=image_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,tc_ZS2_lincomf1f2[[0,1]],tc_ZS2_lincomf1f2[[kk,1]],tc_ZS2_lincomf1f2[[0,(kk+1)]],mu1_lpow,mu2_lpow,mu12_lpow);


    s:=s*l;
  end for;
  "Time for isogeny.", Time(time_isogeny);
  //"finish attack isogeny------------------------";

  lv2tnp_cd:=lv2_img0;
  lv2_x    :=lv2_img_ZS1;
  lv2_y    :=lv2_img_ZS2;

  assert(Is_prod_ell(lv2tnp_to_lv22tnpsq(lv2tnp_cd)));

  //theta transforme to product theta structure.
  lv2tnp_cd,lv2_x,lv2_y:=theta_product_lv2(lv2tnp_cd,lv2_x,lv2_y,zeta_4);

  assert(Is_tehta_prod_ell(lv2tnp_to_lv22tnpsq(lv2tnp_cd)));
  assert(lv2tnp_cd[[0,0]]*lv2tnp_cd[[1,1]] eq lv2tnp_cd[[0,1]]*lv2tnp_cd[[1,0]]);
  assert(lv2_x[[0,0]]*lv2_x[[1,1]] eq lv2_x[[0,1]]*lv2_x[[1,0]]);
  assert(lv2_y[[0,0]]*lv2_y[[1,1]] eq lv2_y[[0,1]]*lv2_y[[1,0]]);
 

  lv2tnp_Ecd_1:=[lv2tnp_cd[[0,0]],lv2tnp_cd[[1,0]]];
  lv2tnp_Ecd_2:=[lv2tnp_cd[[0,0]],lv2tnp_cd[[0,1]]];

  lm_cd_1:=lv2tnp_to_LegendreEll(lv2tnp_Ecd_1);
  lm_cd_2:=lv2tnp_to_LegendreEll(lv2tnp_Ecd_2);

  
  if same_theta_lmd(lm_cd_1,lmd_0) then
    lv2tnp_Ecd:=lv2tnp_Ecd_1;
    lm_cd,sqrt_lm_cd,sqrt_lmm1_cd:=lv2tnp_to_LegendreEll(lv2tnp_Ecd);
    lv2_x_Ecd_1:=[lv2_x[[0,0]],lv2_x[[1,0]]];
    lv2_y_Ecd_1:=[lv2_y[[0,0]],lv2_y[[1,0]]];
    ux_cd:=lv2tc_to_Legendre(lv2tnp_Ecd,sqrt_lm_cd,sqrt_lmm1_cd,lv2_x_Ecd_1);
    uy_cd:=lv2tc_to_Legendre(lv2tnp_Ecd,sqrt_lm_cd,sqrt_lmm1_cd,lv2_y_Ecd_1);
  elif same_theta_lmd(lm_cd_2,lmd_0) then
    lv2tnp_Ecd:=lv2tnp_Ecd_2;
    lm_cd,sqrt_lm_cd,sqrt_lmm1_cd:=lv2tnp_to_LegendreEll(lv2tnp_Ecd);
    lv2_x_Ecd_1:=[lv2_x[[0,0]],lv2_x[[0,1]]];
    lv2_y_Ecd_1:=[lv2_y[[0,0]],lv2_y[[0,1]]];
    ux_cd:=lv2tc_to_Legendre(lv2tnp_Ecd,sqrt_lm_cd,sqrt_lmm1_cd,lv2_x_Ecd_1);
    uy_cd:=lv2tc_to_Legendre(lv2tnp_Ecd,sqrt_lm_cd,sqrt_lmm1_cd,lv2_y_Ecd_1);
  else
    assert(false);
  end if;

  _<x>:=PolynomialRing(GF(p^4));
  Ecd_1:=EllipticCurve(x*(x-1)*(x-lm_cd));

  vx_cd:=Sqrt(ux_cd*(ux_cd-1)*(ux_cd-lm_cd));
  Px:=Ecd_1![ux_cd,vx_cd,1];
  vy_cd:=Sqrt(uy_cd*(uy_cd-1)*(uy_cd-lm_cd));
  Py:=Ecd_1![uy_cd,vy_cd,1];

  assert(IsIsomorphic(Ecd_1,E_0_4));
  _,iso_map:=IsIsomorphic(Ecd_1,E_0_4);

  cand_ker1:=iso_map(Px);
  cand_ker2:=iso_map(Py);
  assert(cand_ker1 in E_0_4);
  assert(cand_ker2 in E_0_4);

  for aut in Automorphisms(E_0_4) do
    ker1:=aut(cand_ker1);
    ker2:=aut(cand_ker2);
    assert(ker1 in E_0_4);
    assert(ker2 in E_0_4);
    if Is_correct_cyclic_isogeny(E_0_4,E_B,N_B,P_A,Q_A,PA_EB,QA_EB,ker1,ker2) then
      _,attacker_kernel:=Is_correct_cyclic_isogeny(E_0_4,E_B,N_B,P_A,Q_A,PA_EB,QA_EB,ker1,ker2);
    end if;
  end for;
  return attacker_kernel;
end function;






//これが最速.
//aux alpha:E_0_4->E_pr.
function MainTorsionAttackKummer_2(E_0_4,E_B,E_pr,N_A,N_B,P_A,Q_A,PA_EB,QA_EB,alpha_P_A,alpha_Q_A,zeta_8)
  "p=      ",p;
  "N_A=    ",N_A;
  "N_B=    ",N_B;
  zeta_4:=zeta_8^2;
  time_total:=Time();
  time_prepare:=Time();
  assert(N_A gt N_B);
  assert(P_A in E_0_4);
  assert(Q_A in E_0_4);
  assert(PA_EB in E_B);
  assert(QA_EB in E_B);
  assert(alpha_P_A in E_pr);
  assert(alpha_Q_A in E_pr);
  lmd_0,lv22tnp_0,lv4tnp_0,E_0_4,j_0,isss_0:=E_to_lmd(E_0_4);
  //"wait for chaing elliptic curve to another isomorphic one."; 
  time_ell_change:=Time();
  lmd_B,lv22tnp_B,lv4tnp_B,E_B2,j_B,isss_B,iso_E_B_2:=E_to_lmd(E_B);
  lmd_pr,lv22tnp_pr,lv4tnp_pr,E_pr,j_pr,isss_pr,iso_E_pr:=E_to_lmd(E_pr);
  //"fin.",Time(time_ell_change);
  assert(isss_0);
  assert(isss_B);
  assert(isss_pr);

  //(N_A,N_A)-isogeny F:E_cd*E_B->E_0*E'_B.
  //basis_KerF:={[alpha(P_A),PA_EB],[alpha(Q_A),QA_EB]}; //in E_cd*E_B.
  //we will call e_1=[alpha(P_A),PA_EB], e_2=[alpha(Q_A),QA_EB].
  //Next we want to calculate  F(0,PA_EB)=(S_1,*), F(0,QA_EB)=(S_2,*), because Ker(phi_B)=<S_1,S_2>.
  S1,S2:=ell_to_torsion_basis_2(E_B2,N_B); //attacker will use.

  _<x>:=PolynomialRing(GF(p^2));
  //E_B_2ex:=EllipticCurve(x*(x-1)*(x-GF(p^2)!lmd_B));
  
  //==========================
  //Theta.
  lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr:=LegendreEll_to_lv2tnp(lmd_pr);
  lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B :=LegendreEll_to_lv2tnp(lmd_B);
  lv2tnp_EprEB:=product_theta(lv2tnp_Epr,lv2tnp_EB);
  lv2tnp_EprEB_2:=prod_lv2tc(E_pr,E_B2,E_pr!0,E_B2!0,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);
  
  f1_Epr:=iso_E_pr(alpha_P_A);    
  f1_EB:=iso_E_B_2(PA_EB); 
  lv2_f1:=prod_lv2tc(E_pr,E_B2,f1_Epr,f1_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);
  
  f2_Epr:=iso_E_pr(alpha_Q_A);    
  f2_EB:=iso_E_B_2(QA_EB);    
  lv2_f2:=prod_lv2tc(E_pr,E_B2,f2_Epr,f2_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);

  f12_Epr:=f1_Epr+f2_Epr;
  f12_EB:=f1_EB+f2_EB;
  lv2_f12:=prod_lv2tc(E_pr,E_B2,f12_Epr,f12_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);

  //------------ZS1.
  ZS1_Epr:=E_pr!0;
  ZS1_EB:=S1;
  lv2_ZS1:=prod_lv2tc(E_pr,E_B2,ZS1_Epr,ZS1_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);

  ZS1f1_Epr:=f1_Epr;
  ZS1f1_EB:=f1_EB+S1;
  lv2_ZS1pf1:=prod_lv2tc(E_pr,E_B2,ZS1f1_Epr,ZS1f1_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);

  ZS1f2_Epr:=f2_Epr;
  ZS1f2_EB:=f2_EB+S1;
  lv2_ZS1pf2:=prod_lv2tc(E_pr,E_B2,ZS1f2_Epr,ZS1f2_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);
  //--------------

  //------------ZS2.
  ZS2_Epr:=E_pr!0;
  ZS2_EB:=S2;
  lv2_ZS2:=prod_lv2tc(E_pr,E_B2,ZS2_Epr,ZS2_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);

  ZS2f1_Epr:=f1_Epr;
  ZS2f1_EB:=f1_EB+S2;
  lv2_ZS2pf1:=prod_lv2tc(E_pr,E_B2,ZS2f1_Epr,ZS2f1_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);

  ZS2f2_Epr:=f2_Epr;
  ZS2f2_EB:=f2_EB+S2;
  lv2_ZS2pf2:=prod_lv2tc(E_pr,E_B2,ZS2f2_Epr,ZS2f2_EB,lv2tnp_Epr,sqrt_lm_pr,sqrt_lmm1_pr,lv2tnp_EB ,sqrt_lm_B ,sqrt_lmm1_B);
  //===========================


  //---------------------------
  lv2_img0      :=lv2tnp_EprEB;
  lv2_img_f1    :=lv2_f1;
  lv2_img_f2    :=lv2_f2;
  lv2_img_f12   :=lv2_f12;
  lv2_img_ZS1   :=lv2_ZS1;
  lv2_img_ZS1pf1:=lv2_ZS1pf1;
  lv2_img_ZS1pf2:=lv2_ZS1pf2;
  lv2_img_ZS2   :=lv2_ZS2;
  lv2_img_ZS2pf1:=lv2_ZS2pf1;
  lv2_img_ZS2pf2:=lv2_ZS2pf2;
  //-----------------------------
 

  "E_0:lmd=",lmd_0;
  "E' :lmd=",lmd_pr;//'
  "E_B:lmd=",lmd_B;
  "";
  "Time for prepare.",Time(time_prepare);
  //=========================================
  time_isogeny:=Time();
  "start attack isogeny-------------------------";
  fac:=decomposition_to_seq(N_A); //the order of composition.
  fac;
  s:=1;

  for i in {1..#fac} do
    "";
    time_thisroop:=Time();
    l:=fac[i];
    kk:=IntegerRing()!(N_A/(s*l));

    //"we calculated degree.", cut_out_first(fac,(i-1));
    //"from now, we calculate degree l=",l;
    "l=",l;
    //"we remain degree.", cut_out_last(fac,(#fac-i));

    assert(s*l*kk eq N_A);
  
    //-----------------------------
    lv2tnp    :=lv2_img0;
    lv2_f1    :=lv2_img_f1;
    lv2_f2    :=lv2_img_f2;
    lv2_f12   :=lv2_img_f12;
    lv2_ZS1   :=lv2_img_ZS1;
    lv2_ZS1pf1:=lv2_img_ZS1pf1;
    lv2_ZS1pf2:=lv2_img_ZS1pf2;
    lv2_ZS2   :=lv2_img_ZS2;
    lv2_ZS2pf1:=lv2_img_ZS2pf1;
    lv2_ZS2pf2:=lv2_img_ZS2pf2;
    //-----------------------------

    //construct kernel of "this" isogeny.
    lv2_e1 :=mult_g2_n2(lv2tnp,kk,lv2_f1);
    lv2_e2 :=mult_g2_n2(lv2tnp,kk,lv2_f2);
    lv2_e12:=mult_g2_n2(lv2tnp,kk,lv2_f12);

    if i ne #fac then  //if not last step.
      lincom_f1f2:=AssociativeArray();
      lincom_f1f2[[(kk+1),0]]:=mult_g2_n2(lv2tnp,kk+1,lv2_f1);
      lincom_f1f2[[0,(kk+1)]]:=mult_g2_n2(lv2tnp,kk+1,lv2_f2);
      lincom_f1f2[[1,kk]]    :=xpke      (lv2tnp,kk  ,lv2_f1,lv2_f2,lv2_f12);
      lincom_f1f2[[kk,1]]    :=xpke      (lv2tnp,kk  ,lv2_f2,lv2_f1,lv2_f12);
    end if;

   
    //ZS1----------------------------------
    //"wait for calculate about ZS1.";
    tc_ZS1_lincomf1f2:=AssociativeArray();
    tc_ZS1_lincomf1f2[[0,0]]:=lv2_ZS1;
    tc_ZS1_lincomf1f2[[kk,0]]:=xpke(lv2tnp,kk,lv2_ZS1,lv2_f1,lv2_ZS1pf1);
    tc_ZS1_lincomf1f2[[0,kk]]:=xpke(lv2tnp,kk,lv2_ZS1,lv2_f2,lv2_ZS1pf2);
    //ZS2----------------------------------
    //"wait for calculate about ZS2.";
    tc_ZS2_lincomf1f2:=AssociativeArray();
    tc_ZS2_lincomf1f2[[0,0]]:=lv2_ZS2;
    tc_ZS2_lincomf1f2[[kk,0]]:=xpke(lv2tnp,kk,lv2_ZS2,lv2_f1,lv2_ZS2pf1);
    tc_ZS2_lincomf1f2[[0,kk]]:=xpke(lv2tnp,kk,lv2_ZS2,lv2_f2,lv2_ZS2pf2);
    //================================

    //img of 0.
    lv2_img0,mu1_lpow,mu2_lpow,mu12_lpow:=codomain_lv2tnp(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l);
    //---------
    if i ne #fac then  //if not last step.
      //img of f_1.
      lv2_img_f1:=image_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,lv2_f1,lincom_f1f2[[(kk+1),0]],lincom_f1f2[[1,kk]],mu1_lpow,mu2_lpow,mu12_lpow);
      //img of f_2.
      lv2_img_f2:=image_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,lv2_f2,lincom_f1f2[[kk,1]],lincom_f1f2[[0,(kk+1)]],mu1_lpow,mu2_lpow,mu12_lpow);
      //img of f_12.
      lv2_img_f12:=Add_lv2(lv2_img0,lv2_img_f1,lv2_img_f2);
    end if;
    //--------------
    //img of ZS1.
    lv2_img_ZS1:=image_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,tc_ZS1_lincomf1f2[[0,0]],tc_ZS1_lincomf1f2[[kk,0]],tc_ZS1_lincomf1f2[[0,kk]],mu1_lpow,mu2_lpow,mu12_lpow);

    if i ne #fac then //if not last step.
      //img of ZS1+f_1.
      lv2_img_ZS1pf1:=Add_lv2(lv2_img0,lv2_img_ZS1,lv2_img_f1);
      //img of ZS1+f_2.
      lv2_img_ZS1pf2:=CompatibleAdd(lv2_img0,lv2_img_ZS1,lv2_img_f2,lv2_img_ZS1pf1,lv2_img_f12);
    end if;

    //--------------
    //img of ZS2.
    lv2_img_ZS2:=image_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,tc_ZS2_lincomf1f2[[0,0]],tc_ZS2_lincomf1f2[[kk,0]],tc_ZS2_lincomf1f2[[0,kk]],mu1_lpow,mu2_lpow,mu12_lpow);

    if i ne #fac then //if not last step.
      //img of ZS2+f_1.
      lv2_img_ZS2pf1:=Add_lv2(lv2_img0,lv2_img_ZS2,lv2_img_f1);
      //img of ZS2+f_2.
      lv2_img_ZS2pf2:=CompatibleAdd(lv2_img0,lv2_img_ZS2,lv2_img_f2,lv2_img_ZS2pf1,lv2_img_f12);
    end if;

    "time_roop.",Time(time_thisroop);

    s:=s*l;
  end for;
  "Time for isogeny.", Time(time_isogeny);
  //"finish attack isogeny------------------------";

  lv2tnp_cd:=lv2_img0;
  lv2_x    :=lv2_img_ZS1;
  lv2_y    :=lv2_img_ZS2;

  assert(Is_prod_ell(lv2tnp_to_lv22tnpsq(lv2tnp_cd)));

  //theta transforme to product theta structure.
  lv2tnp_cd,lv2_x,lv2_y:=theta_product_lv2(lv2tnp_cd,lv2_x,lv2_y,zeta_4);

  assert(Is_tehta_prod_ell(lv2tnp_to_lv22tnpsq(lv2tnp_cd)));
  assert(lv2tnp_cd[[0,0]]*lv2tnp_cd[[1,1]] eq lv2tnp_cd[[0,1]]*lv2tnp_cd[[1,0]]);
  assert(lv2_x[[0,0]]*lv2_x[[1,1]] eq lv2_x[[0,1]]*lv2_x[[1,0]]);
  assert(lv2_y[[0,0]]*lv2_y[[1,1]] eq lv2_y[[0,1]]*lv2_y[[1,0]]);
 

  lv2tnp_Ecd_1:=[lv2tnp_cd[[0,0]],lv2tnp_cd[[1,0]]];
  lv2tnp_Ecd_2:=[lv2tnp_cd[[0,0]],lv2tnp_cd[[0,1]]];

  lm_cd_1:=lv2tnp_to_LegendreEll(lv2tnp_Ecd_1);
  lm_cd_2:=lv2tnp_to_LegendreEll(lv2tnp_Ecd_2);

  
  if same_theta_lmd(lm_cd_1,lmd_0) then
    lv2tnp_Ecd:=lv2tnp_Ecd_1;
    lm_cd,sqrt_lm_cd,sqrt_lmm1_cd:=lv2tnp_to_LegendreEll(lv2tnp_Ecd);
    lv2_x_Ecd_1:=[lv2_x[[0,0]],lv2_x[[1,0]]];
    lv2_y_Ecd_1:=[lv2_y[[0,0]],lv2_y[[1,0]]];
    ux_cd:=lv2tc_to_Legendre(lv2tnp_Ecd,sqrt_lm_cd,sqrt_lmm1_cd,lv2_x_Ecd_1);
    uy_cd:=lv2tc_to_Legendre(lv2tnp_Ecd,sqrt_lm_cd,sqrt_lmm1_cd,lv2_y_Ecd_1);
  elif same_theta_lmd(lm_cd_2,lmd_0) then
    lv2tnp_Ecd:=lv2tnp_Ecd_2;
    lm_cd,sqrt_lm_cd,sqrt_lmm1_cd:=lv2tnp_to_LegendreEll(lv2tnp_Ecd);
    lv2_x_Ecd_1:=[lv2_x[[0,0]],lv2_x[[0,1]]];
    lv2_y_Ecd_1:=[lv2_y[[0,0]],lv2_y[[0,1]]];
    ux_cd:=lv2tc_to_Legendre(lv2tnp_Ecd,sqrt_lm_cd,sqrt_lmm1_cd,lv2_x_Ecd_1);
    uy_cd:=lv2tc_to_Legendre(lv2tnp_Ecd,sqrt_lm_cd,sqrt_lmm1_cd,lv2_y_Ecd_1);
  else
    assert(false);
  end if;

  _<x>:=PolynomialRing(GF(p^4));
  Ecd_1:=EllipticCurve(x*(x-1)*(x-lm_cd));

  vx_cd:=Sqrt(ux_cd*(ux_cd-1)*(ux_cd-lm_cd));
  Px:=Ecd_1![ux_cd,vx_cd,1];
  vy_cd:=Sqrt(uy_cd*(uy_cd-1)*(uy_cd-lm_cd));
  Py:=Ecd_1![uy_cd,vy_cd,1];

  assert(IsIsomorphic(Ecd_1,E_0_4));
  _,iso_map:=IsIsomorphic(Ecd_1,E_0_4);

  cand_ker1:=iso_map(Px);
  cand_ker2:=iso_map(Py);
  assert(cand_ker1 in E_0_4);
  assert(cand_ker2 in E_0_4);

  for aut in Automorphisms(E_0_4) do
    ker1:=aut(cand_ker1);
    ker2:=aut(cand_ker2);
    assert(ker1 in E_0_4);
    assert(ker2 in E_0_4);
    if Is_correct_cyclic_isogeny(E_0_4,E_B,N_B,P_A,Q_A,PA_EB,QA_EB,ker1,ker2) then
      _,attacker_kernel:=Is_correct_cyclic_isogeny(E_0_4,E_B,N_B,P_A,Q_A,PA_EB,QA_EB,ker1,ker2);
    end if;
  end for;

  "Time total.", Time(time_total);
  return attacker_kernel;
end function;




