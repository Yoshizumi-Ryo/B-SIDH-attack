//level 2.

//==========================-
//Not functions.

lv2keys:={[a,b]:a,b in {0,1}};


Riemann_Relation_list:={};
for comp in CartesianPower(lv2keys,5) do
  chi:=comp[1];
  ii:=comp[2];
  jj:=comp[3];
  kk:=comp[4];
  ll:=[(ii[1]+jj[1]+kk[1])mod 2,(ii[2]+jj[2]+kk[2])mod 2];
  mm:=comp[5];
  id:=[(mm[1]-ii[1])mod 2,(mm[2]-ii[2])mod 2];
  jd:=[(mm[1]-jj[1])mod 2,(mm[2]-jj[2])mod 2];
  kd:=[(mm[1]-kk[1])mod 2,(mm[2]-kk[2])mod 2];
  ld:=[(mm[1]-ll[1])mod 2,(mm[2]-ll[2])mod 2];
  Riemann_Relation_list join:={[chi,ii,jj,kk,ll,id,jd,kd,ld]};
end for;




represent_set:=
 {[[1,0],[1,1],[1,0],[1,0],[1,1],[0,1],[0,0],[0,0],[0,1]],[[0,1],[0,1],[1,1],[1,1],[0,1],[1,0],[0,0],[0,0],[1,0]],[[1,1],[1,0],[1,0],[0,0],[0,0],[0,1],[0,1],[1,1],[1,1]],[[0,0],[1,0],[0,0],[1,1],[0,1],[1,0],[0,0],[1,1],[0,1]],[[1,1],[1,0],[1,0],[1,0],[1,0],[1,1],[1,1],[1,1],[1,1]],[[1,0],[1,0],[1,0],[0,0],[0,0],[0,1],[0,1],[1,1],[1,1]],[[0,1],[0,1],[0,1],[1,1],[1,1],[1,0],[1,0],[0,0],[0,0]],[[1,0],[0,0],[0,1],[1,0],[1,1],[0,0],[0,1],[1,0],[1,1]],[[0,0],[1,1],[1,1],[1,1],[1,1],[1,0],[1,0],[1,0],[1,0]],[[1,0],[0,0],[1,0],[0,1],[1,1],[0,1],[1,1],[0,0],[1,0]],[[0,1],[0,1],[0,1],[0,0],[0,0],[1,0],[1,0],[1,1],[1,1]],[[1,0],[1,1],[1,1],[1,1],[1,1],[0,0],[0,0],[0,0],[0,0]],[[0,1],[1,0],[0,0],[0,1],[1,1],[0,0],[1,0],[1,1],[0,1]],[[0,0],[0,1],[0,0],[1,0],[1,1],[1,1],[1,0],[0,0],[0,1]],[[0,0],[1,1],[0,0],[1,1],[0,0],[1,0],[0,1],[1,0],[0,1]]
};




//============================
//lv2 theta on elliptic curve.



function eq_lv2tc(lv2tc_1,lv2tc_2)
  ratio:=lv2tc_1[[0,0]]/lv2tc_2[[0,0]];
  ct:=0;
  for key in lv2keys do
    if lv2tc_1[key] eq ratio*lv2tc_2[key] then
      ct+:=1;
    end if;
  end for;
  if ct eq 4 then
    return true;
  else
    return false;
  end if;
end function;


//from y^2=x(x-1)(x-lm) to level 2 theta null point.
function LegendreEll_to_lv2tnp(lm)
  _<x>:=PolynomialRing(Parent(lm));
  sqrt_lm:=RootsInSplittingField(x^2-lm)[1][1];
  sqrt_lmm1:=RootsInSplittingField(x^2-(lm-1))[1][1];
  thnp0_sq:=sqrt_lm;
  thnp1_sq:=sqrt_lmm1;
  thnp2_sq:=1;
  thnp3_sq:=0;
  lv2tnp:=AssociativeArray();
  lv2tnp0:=thnp0_sq+thnp2_sq;
  lv2tnp1:=thnp1_sq;
  lv2tnp:=[lv2tnp0,lv2tnp1];
  return lv2tnp,sqrt_lm,sqrt_lmm1;
end function;


//the inverse function.
function lv2tnp_to_LegendreEll(lv2tnp)
  lm:=((lv2tnp[1]^2+lv2tnp[2]^2)/(lv2tnp[1]^2-lv2tnp[2]^2))^2;
  r:=lv2tnp[1]/lv2tnp[2];
  sqrt_lm:=((lm-1)*r^2-lm-1)/2;
  sqrt_lmm1:=(sqrt_lm+1)/r;
  assert(sqrt_lm^2 eq lm);
  assert(sqrt_lmm1^2 eq lm-1);
  return lm,sqrt_lm,sqrt_lmm1;
end function;


//for P in E, we give level 2 theta coordinate.
function Legendre_to_lv2tc(lv2tnp,sqrt_lm,sqrt_lmm1,P)
  if P eq Parent(P)!0 then
    return lv2tnp;
  end if;
  u:=P[1];
  lm:=lv2tnp_to_LegendreEll(lv2tnp);
  //assert(v^2 eq u*(u-1)*(u-lm));
  assert(sqrt_lm^2 eq lm);
  assert(sqrt_lmm1^2 eq lm-1);
  thc0_sq:=sqrt_lm*(u-1);
  thc1_sq:=sqrt_lmm1*u;
  thc2_sq:=sqrt_lm*thc0_sq-sqrt_lmm1*thc1_sq;
  thc3_sq:=sqrt_lm*thc1_sq-sqrt_lmm1*thc0_sq;
  return [1,(thc1_sq+thc3_sq)/(thc0_sq+thc2_sq)];
end function;


//the inverse function.
function lv2tc_to_Legendre(lv2tnp,sqrt_lm,sqrt_lmm1,lv2tc)
  lm:=sqrt_lm^2;
  r:=lv2tc[2]/lv2tc[1];
  u:=(sqrt_lm*sqrt_lmm1+r*(lm+sqrt_lm))/(r*(1+sqrt_lm)-sqrt_lmm1);
  return u;
end function;


//from lv2 theta null point to squared lv(2,2) theta null point.
function lv2tnp_to_lv22tnpsq(lv2tnp)
  lv22tnpsq:=AssociativeArray();
  for key in lv22keys do
    lv22tnpsq[key]:=0;
    for beta in {[0,0],[0,1],[1,0],[1,1]} do
      lv22tnpsq[key]+:=(-1)^(key[1]*beta[1]+key[2]*beta[2])*lv2tnp[[(key[3]+beta[1])mod 2,(key[4]+beta[2])mod 2]]*lv2tnp[beta];
    end for;
  end for;
  return lv22tnpsq;
end function;



//level 2 theta null point of elliptic product.
function product_theta(lv2tc_1,lv2tc_2)
  lv2tc:=AssociativeArray();
  lv2tc[[0,0]]:=lv2tc_1[1]*lv2tc_2[1];
  lv2tc[[0,1]]:=lv2tc_1[1]*lv2tc_2[2];
  lv2tc[[1,0]]:=lv2tc_1[2]*lv2tc_2[1];
  lv2tc[[1,1]]:=lv2tc_1[2]*lv2tc_2[2];
  return lv2tc;
end function;


 

//about addition.==================================


//product term of Riemann relation.
function DiffAdd_product_term(chi,lv2_x,lv2_y,i,j)
  sum:=0;
  for t in {[a,b]:a,b in {0,1}} do 
    sign:=(-1)^(chi[1]*t[1]+chi[2]*t[2]);
    ipt:=[(i[1]+t[1]) mod 2,(i[2]+t[2]) mod 2];
    jpt:=[(j[1]+t[1]) mod 2,(j[2]+t[2]) mod 2];
    sum+:=sign*lv2_x[ipt]*lv2_y[jpt];
  end for;
  return sum;
end function;




//normal addition.
function Add_lv2(lv2tnp,lv2_x,lv2_y)
  //assert(not(Is_prod_ell(lv2tnp_to_lv22tnpsq(lv2tnp))));
  Kxx<r,X01,X10,X11,Y01,Y10,Y11>:=PolynomialRing(GF(p^4),7);
  
  X:=AssociativeArray();
  Y:=AssociativeArray();
  X[[0,0]] := GF(p^4)!1;
  X[[0,1]] := X01;
  X[[1,0]] := X10;
  X[[1,1]] := X11;
  Y[[0,0]] := GF(p^4)!1;
  Y[[0,1]] := Y01;
  Y[[1,0]] := Y10;
  Y[[1,1]] := Y11;
  
  f:=AssociativeArray();
  for chirr in Riemann_Relation_list do
    chi:=chirr[1];
    i:=chirr[2];
    j:=chirr[3];
    k:=chirr[4];
    l:=chirr[5];
    id:=chirr[6];
    jd:=chirr[7];
    kd:=chirr[8];
    ld:=chirr[9];
    LHSfst:=&+[r*((-1)^(chi[1]*t[1]+chi[2]*t[2]))*(X[[(i[1]+t[1])mod 2,(i[2]+t[2])mod 2]])*(Y[[(j[1]+t[1])mod 2,(j[2]+t[2])mod 2]]):t in lv2keys];
    f[chirr]:=LHSfst*Kxx!DiffAdd_product_term(chi,lv2tnp,lv2tnp,k,l)-Kxx!DiffAdd_product_term(chi,lv2_x,lv2_x,id,jd)*Kxx!DiffAdd_product_term(chi,lv2_y,lv2_y,kd,ld);
  end for;
  polysis:= {f[key]:key in Keys(f)};
  //#Keys(f);
  //#Riemann_Relation_list;
  //#polysis;
  I:=ideal<Kxx|polysis>;
  Groebner(I);
  V:=Variety(I);
  assert(#V ne 0);
  sol:=V[1];
  lv2_xpy:=AssociativeArray();
  lv2_xpy[[0,0]]:=1;
  lv2_xpy[[0,1]]:=sol[2];
  lv2_xpy[[1,0]]:=sol[3];
  lv2_xpy[[1,1]]:=sol[4];
  lv2_xmy:=AssociativeArray();
  lv2_xmy[[0,0]]:=1;
  lv2_xmy[[0,1]]:=sol[5];
  lv2_xmy[[1,0]]:=sol[6];
  lv2_xmy[[1,1]]:=sol[7];
  return lv2_xpy,lv2_xmy;
end function;


//cf.[LR12].section3.2.1.
function CompatibleAdd(lv2tnp,lv2_x,lv2_z,lv2_xpy,lv2_ypz)
  lv2_xpz_1,lv2_xpz_2:=Add_lv2(lv2tnp,lv2_x,lv2_z);
  lv2tc_1,lv2tc_2:=Add_lv2(lv2tnp,lv2_xpy,lv2_xpz_1);
  lv2tc_3,lv2tc_4:=Add_lv2(lv2tnp,lv2_xpy,lv2_xpz_2);
  if eq_lv2tc(lv2tc_1,lv2_ypz) or eq_lv2tc(lv2tc_2,lv2_ypz) then
    return lv2_xpz_2;
  elif eq_lv2tc(lv2tc_3,lv2_ypz) or eq_lv2tc(lv2tc_4,lv2_ypz) then
    return lv2_xpz_1;
  else
    assert(false);
  end if;
end function;


      
//正しいが遅い. 
function DiffAdd_lv2_linsys(lv2tnp,lv2_x,lv2_y,lv2_xmy)
  Kxx<X00,X01,X10,X11>:=PolynomialRing(GF(p^4),4);
  X:=AssociativeArray();
  X[[0,0]] := X00;
  X[[0,1]] := X01;
  X[[1,0]] := X10;
  X[[1,1]] := X11;
  f:=AssociativeArray();
  for chirr in Riemann_Relation_list do
    chi:=chirr[1];
    i:=chirr[2];
    j:=chirr[3];
    k:=chirr[4];
    l:=chirr[5];
    id:=chirr[6];
    jd:=chirr[7];
    kd:=chirr[8];
    ld:=chirr[9];
    LHSfst:=&+[((-1)^(chi[1]*t[1]+chi[2]*t[2]))*(X[[(i[1]+t[1])mod 2,(i[2]+t[2])mod 2]])*(lv2_xmy[[(j[1]+t[1])mod 2,(j[2]+t[2])mod 2]]):t in lv2keys];
    f[chirr]:=LHSfst*Kxx!DiffAdd_product_term(chi,lv2tnp,lv2tnp,k,l)-Kxx!DiffAdd_product_term(chi,lv2_x,lv2_x,id,jd)*Kxx!DiffAdd_product_term(chi,lv2_y,lv2_y,kd,ld);
  end for;
  polysis:= {f[key]:key in Keys(f)};
  //#polysis;
  I:=ideal<Kxx|polysis>;
  Groebner(I);
  V:=Variety(I);
  //"#V",#V;
  assert(#V eq 1);
  lv2_xpy:=AssociativeArray();
  lv2_xpy[[0,0]]:=V[1][1];
  lv2_xpy[[0,1]]:=V[1][2];
  lv2_xpy[[1,0]]:=V[1][3];
  lv2_xpy[[1,1]]:=V[1][4];
  return lv2_xpy;
end function;



//現時点, 最速のdifferential addition of g=2,n=2.
function DiffAdd_lv2_linsys_2(lv2tnp,lv2_x,lv2_y,lv2_xmy)
  Kxx<X00,X01,X10,X11>:=PolynomialRing(GF(p^4),4);
  X:=AssociativeArray();
  X[[0,0]] := X00;
  X[[0,1]] := X01;
  X[[1,0]] := X10;
  X[[1,1]] := X11;
  polysis:={};
  for chirr in represent_set do
    chi:=chirr[1];
    ii:=chirr[2];
    jj:=chirr[3];
    kk:=chirr[4];
    ll:=chirr[5];
    id:=chirr[6];
    jd:=chirr[7];
    kd:=chirr[8];
    ld:=chirr[9];
    LHSfst:=&+[((-1)^(chi[1]*t[1]+chi[2]*t[2]))*(X[[(ii[1]+t[1])mod 2,(ii[2]+t[2])mod 2]])*(lv2_xmy[[(jj[1]+t[1])mod 2,(jj[2]+t[2])mod 2]]):t in lv2keys];
    f:=LHSfst*Kxx!DiffAdd_product_term(chi,lv2tnp,lv2tnp,kk,ll)-Kxx!DiffAdd_product_term(chi,lv2_x,lv2_x,id,jd)*Kxx!DiffAdd_product_term(chi,lv2_y,lv2_y,kd,ld);
    polysis join:={f};
  end for;
  assert(#polysis eq 15);
  I:=ideal<Kxx|polysis>;
  Groebner(I);
  V:=Variety(I);
  assert(#V eq 1);
  lv2_xpy:=AssociativeArray();
  lv2_xpy[[0,0]]:=V[1][1];
  lv2_xpy[[0,1]]:=V[1][2];
  lv2_xpy[[1,0]]:=V[1][3];
  lv2_xpy[[1,1]]:=V[1][4];
  return lv2_xpy;
end function;




function Double_g2n2(lv2tnp,lv2_x)
  return DiffAdd_lv2_linsys_2(lv2tnp,lv2_x,lv2_x,lv2tnp);
end function;



//return x+y+z.
function Extended_Add(lv2tnp,lv2_x,lv2_y,lv2_z,lv2_xpy,lv2_ypz,lv2_xpz)
  Kxx<X00,X01,X10,X11>:=PolynomialRing(GF(p^4),4);
  X:=AssociativeArray();
  X[[0,0]] := X00;
  X[[0,1]] := X01;
  X[[1,0]] := X10;
  X[[1,1]] := X11;
  polysis:={};
  for chirr in represent_set do
    chi:=chirr[1];
    ii:=chirr[2];
    jj:=chirr[3];
    kk:=chirr[4];
    ll:=chirr[5];
    id:=chirr[6];
    jd:=chirr[7];
    kd:=chirr[8];
    ld:=chirr[9];
    LHSfst:=&+[((-1)^(chi[1]*t[1]+chi[2]*t[2]))*(X[[(ii[1]+t[1])mod 2,(ii[2]+t[2])mod 2]])*(lv2_x[[(jj[1]+t[1])mod 2,(jj[2]+t[2])mod 2]]):t in lv2keys];
    f:=LHSfst*Kxx!DiffAdd_product_term(chi,lv2_y,lv2_z,kk,ll)-Kxx!DiffAdd_product_term(chi,lv2tnp,lv2_ypz,id,jd)*Kxx!DiffAdd_product_term(chi,lv2_xpz,lv2_xpy,kd,ld);
    polysis join:={f};
  end for;
  assert(#polysis eq 15);
  I:=ideal<Kxx|polysis>;
  Groebner(I);
  V:=Variety(I);
  //#V;
  assert(#V eq 1);
  lv2_xpypz:=AssociativeArray();
  lv2_xpypz[[0,0]]:=V[1][1];
  lv2_xpypz[[0,1]]:=V[1][2];
  lv2_xpypz[[1,0]]:=V[1][3];
  lv2_xpypz[[1,1]]:=V[1][4];
  return lv2_xpypz;
end function;
  



//Ladder.
function Ladder_g2n2(lv2tnp,k,lv2_x)
  bit_k:=IntegerToSequence(k,2);
  x:=lv2_x;
  y:=DiffAdd_lv2_linsys_2(lv2tnp,lv2_x,lv2_x,lv2tnp);
  xmy:=lv2_x;
  for i in [1..(#bit_k-1)] do
    bit_i:=bit_k[#bit_k-i];  //上からj番目
    if bit_i eq 1 then
      x_0:=DiffAdd_lv2_linsys_2(lv2tnp,x,y,xmy);
      y_0:=DiffAdd_lv2_linsys_2(lv2tnp,y,y,lv2tnp);
      x:=x_0;
      y:=y_0;
    end if;
    if bit_i eq 0 then
      x_0:=DiffAdd_lv2_linsys_2(lv2tnp,x,x,lv2tnp);
      y_0:=DiffAdd_lv2_linsys_2(lv2tnp,x,y,xmy);
      x:=x_0;
      y:=y_0;
    end if;
  end for;
  return x;
end function;



//k倍の計算.現時点で最速.
function mult_g2_n2(lv2tnp,k,lv2_x)
  assert(k ge 0);
  if k eq 0 then
    return lv2tnp;
  elif k eq 1 then
    return lv2_x;  
  elif ((k ge 2) and (k le 5)) then
    lv2_km2x:=lv2tnp;
    lv2_km1x:=lv2_x;
    for i in [2..k] do
      lv2_kx:=DiffAdd_lv2_linsys_2(lv2tnp,lv2_km1x,lv2_x,lv2_km2x);
      lv2_km2x:=lv2_km1x;
      lv2_km1x:=lv2_kx;
    end for;
    return lv2_kx;
  else
    return Ladder_g2n2(lv2tnp,k,lv2_x);  
  end if;
end function;




//compute kx+y from x,y,x-y on lv2.
function ThreePtLadder_g2n2(alnp,k,alx,aly,alxmy)
  bit_k:=IntegerToSequence(k,2);
  X:=alx;
  Y:=aly;
  Z:=alxmy;
  for i in [1..(#bit_k)] do
    if bit_k[i] eq 1 then
      X0:=Double_g2n2(alnp,X); 
      Y:=DiffAdd_lv2_linsys_2(alnp,X,Y,Z);
      X:=X0;
    else
      X0:=Double_g2n2(alnp,X);
      Z:=DiffAdd_lv2_linsys_2(alnp,X,Z,Y);
      X:=X0;
    end if;
  end for;
  return Y;
end function;




//現時点最速.compute x+ke.
function xpke(lv2tnp,k,lv2_x,lv2_e,lv2_xpe)
  if k eq 0 then
    return lv2_x; //k=0.
  elif k eq 1 then
    return lv2_xpe;  //k=1.
  elif k ge 2 then
    return ThreePtLadder_g2n2(lv2tnp,k-1,lv2_e,lv2_xpe,lv2_x);  
  end if;
end function;




//x+ke. 愚直な計算.
function xpke_2(lv2tnp,lv2_x,k,lv2_e,lv2_xpe)
  if k eq 0 then
    return lv2_x;
  elif k eq 1 then
    return lv2_xpe;
  elif k ge 2 then
    lv2_xpkm2e:=lv2_x;  //x+(k-2)e
    lv2_xpkm1e:=lv2_xpe;  //x+(k-1)e
    for i in [2..k] do
      lv2_xpke:=DiffAdd_lv2_linsys_2(lv2tnp,lv2_xpkm1e,lv2_e,lv2_xpkm2e); //x+ke
      lv2_xpkm2e:=lv2_xpkm1e;
      lv2_xpkm1e:=lv2_xpke;
    end for;
  end if;
  return lv2_xpke;
end function;




//for codomain theta null point.======================

//calculate excellent  l-torsion point e on g=2,n=2.
function calculate_mu_lpow(lv2tnp,lv2_e,l)
  _,l_d:=IsDivisibleBy(l-1,2);  //l=2l'+1.
  for i_0 in lv2keys do
    denom:=mult_g2_n2(lv2tnp,l_d+1,lv2_e)[i_0];
    if denom ne 0 then
      lm_lpow:=mult_g2_n2(lv2tnp,l_d,lv2_e)[i_0]/denom;
      return lm_lpow;
    end if;
  end for;
end function;


//calculate linear combination k1e1+k2e2 for 0<=k1,k2<=(l-1).
function lincom_e1e2_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l)
  lincom:=AssociativeArray();
  lincom[[0,0]]:=lv2tnp;
  lincom[[0,1]]:=lv2_e2;
  for k2 in [2..(l-1)] do
    lincom[[0,k2]]:=DiffAdd_lv2_linsys_2(lv2tnp,lincom[[0,k2-1]],lv2_e2,lincom[[0,k2-2]]);
  end for;
  lincom[[1,0]]:=lv2_e1;
  lincom[[1,1]]:=lv2_e12;
  for k2 in [2..(l-1)] do
    lincom[[1,k2]]:=DiffAdd_lv2_linsys_2(lv2tnp,lincom[[1,k2-1]],lv2_e2,lincom[[1,k2-2]]);
  end for;
  for k2 in [0..(l-1)] do
    for k1 in [2..(l-1)] do
      lincom[[k1,k2]]:=DiffAdd_lv2_linsys_2(lv2tnp,lincom[[k1-1,k2]],lv2_e1,lincom[[k1-2,k2]]);
    end for;
  end for;
  assert(#Keys(lincom) eq l^2);
  return lincom;
end function;



function codomain_lv2tnp(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l)

  time_lincom:=Time();
  lincom:=lincom_e1e2_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l);
  "time_limcom.  ",Time(time_lincom);

  mu1_lpow:=calculate_mu_lpow(lv2tnp,lv2_e1,l);
  mu2_lpow:=calculate_mu_lpow(lv2tnp,lv2_e2,l);
  mu12_lpow:=calculate_mu_lpow(lv2tnp,lv2_e12,l);
  lv2tnp_cd:=AssociativeArray();
  for j in lv2keys do
    lv2tnp_cd[j]:=0;
    for k1 in [0..(l-1)] do
      for k2 in [0..(l-1)] do
        lv2tnp_cd[j]+:=mu1_lpow^(k1^2-k1*k2)*mu2_lpow^(k2^2-k1*k2)*mu12_lpow^(k1*k2)*lincom[[k1,k2]][j]^l;
      end for;
    end for;
  end for;
  return lv2tnp_cd,mu1_lpow,mu2_lpow,mu12_lpow;
end function;


//for iamge of point.===============================


function calculate_nu_lpow(lv2tnp,lv2_x,lv2_e,lv2_xpe,l,mu_lpow)
  lv2_xple:=xpke(lv2tnp,l,lv2_x,lv2_e,lv2_xpe);
  nu_lpow:=lv2_x[[0,0]]/((mu_lpow^(l-1))*lv2_xple[[0,0]]);
  return nu_lpow;
end function;



function xplincom_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,lv2_x,lv2_xpe1,lv2_xpe2)
  xplincom:=AssociativeArray();
  xplincom[[0,0]]:=lv2_x;
  xplincom[[1,0]]:=lv2_xpe1;
  xplincom[[0,1]]:=lv2_xpe2;
  xplincom[[1,1]]:=Extended_Add(lv2tnp,lv2_x,lv2_e1,lv2_e2,lv2_xpe1,lv2_e12,lv2_xpe2);
  for k2 in [2..(l-1)] do
    xplincom[[0,k2]]:=DiffAdd_lv2_linsys_2(lv2tnp,xplincom[[0,k2-1]],lv2_e2,xplincom[[0,k2-2]]);
  end for;
  for k2 in [2..(l-1)] do
    xplincom[[1,k2]]:=DiffAdd_lv2_linsys_2(lv2tnp,xplincom[[1,k2-1]],lv2_e2,xplincom[[1,k2-2]]);
  end for;
  for k2 in [0..(l-1)] do
    for k1 in [2..(l-1)] do
      xplincom[[k1,k2]]:=DiffAdd_lv2_linsys_2(lv2tnp,xplincom[[k1-1,k2]],lv2_e1,xplincom[[k1-2,k2]]);
    end for;
  end for;
  assert(#Keys(xplincom) eq l^2);
  return xplincom;
end function;

  

function image_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,lv2_x,lv2_xpe1,lv2_xpe2,mu1_lpow,mu2_lpow,mu12_lpow)
  nu1_lpow:=calculate_nu_lpow(lv2tnp,lv2_x,lv2_e1,lv2_xpe1,l,mu1_lpow);
  nu2_lpow:=calculate_nu_lpow(lv2tnp,lv2_x,lv2_e2,lv2_xpe2,l,mu2_lpow);

  time_xplincom:=Time();
  xpKer:=xplincom_g2n2(lv2tnp,lv2_e1,lv2_e2,lv2_e12,l,lv2_x,lv2_xpe1,lv2_xpe2);
  "time_xplincom.",Time(time_xplincom);

  lv2_fx:=AssociativeArray();
  for j in lv2keys do
    lv2_fx[j]:=0;
    for k1 in {0..(l-1)} do
      for k2 in {0..(l-1)} do
        excell_coff:=mu1_lpow^(k1*(k1-k2-1))*mu2_lpow^(k2*(k2-k1-1))*mu12_lpow^(k1*k2)*nu1_lpow^k1*nu2_lpow^k2;
        lv2_fx[j]+:=excell_coff*((xpKer[[k1,k2]][j])^l);
      end for;
    end for;
  end for;
  return lv2_fx;
end function;




//other useful functions.=============================

procedure get_order_lv2(lv2tnp,lv2_x)
  for ct in {1..100} do
    if eq_lv2tc(lv2tnp,mult_g2_n2(lv2tnp,ct,lv2_x)) then
      ct;
      break ct;
    end if;
  end for;
end procedure;




function IsOrder_g2lv2(lv2tnp,k,lv2_x)
  lv2_kx:=mult_g2_n2(lv2tnp,k,lv2_x);
  return eq_lv2tc(lv2tnp,lv2_kx);
end function;




//====================================