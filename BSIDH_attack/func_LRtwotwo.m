//[LR22]algotirhm.



//Here, lv4tc is l-torsion point.
function excellent_torsion(lv4tnp,lv4_e,l)
  _,l_d:=IsDivisibleBy(l-1,2);  //l=2l'-1.
  inv_lv4_e:=inverse_element(lv4_e);
  for i_0 in lv4keys do
    if (opt_mult(lv4tnp,l_d+1,lv4_e)[i_0] ne 0) then
      lm_lpow:=opt_mult(lv4tnp,l_d,inv_lv4_e)[i_0]/opt_mult(lv4tnp,l_d+1,lv4_e)[i_0];
      break i_0;
    end if;
  end for;
  _<x>:=PolynomialRing(Parent(lm_lpow));
  lm:=RootsInSplittingField(x^l-lm_lpow)[1][1];
  for i in lv4keys do
    lv4_e[i]*:=lm;
  end for;
  return lv4_e;
end function;


//for l-torsion point e.
function calculate_mu_lpow(lv4tnp,lv4_e,l)
  _,l_d:=IsDivisibleBy(l-1,2);  //l=2l'+1.
  inv_lv4_e:=inverse_element(lv4_e);
  for i_0 in lv4keys do
    denom:=opt_mult(lv4tnp,l_d+1,lv4_e)[i_0];
    if denom ne 0 then
      lm_lpow:=opt_mult(lv4tnp,l_d,inv_lv4_e)[i_0]/denom;
      return lm_lpow;
    end if;
  end for;
end function;









//see [LR22]Def3.7
function check_excellent(lv4tnp,lv4tc,l)
  _,l_d:=IsDivisibleBy(l-1,2); 
  assert(l eq (2*l_d+1));
  inv_lv4tc:=inverse_element(lv4tc);
  return eq_Assoc(mult3(lv4tnp,l_d+1,lv4tc),mult3(lv4tnp,l_d,inv_lv4tc));
end function;




function check_excellent2(lv4tnp,lv4tc,l)
  return eq_Assoc(mult3(lv4tnp,l,lv4tc),lv4tnp);
end function;






//P in A[l].
//x is any point in A.
//[LR22]Alg1.
function excellent_pt(lv4tnp,lv4_P,lv4_x,lv4_Ppx,l)
  lv4_P:=excellent_torsion(lv4tnp,lv4_P,l);
  assert(check_excellent(lv4tnp,lv4_P,l));
  lv4_lPpx:=ThreePtLadder_plus(lv4tnp,l,lv4_P,lv4_x,lv4_Ppx); //lP+x.
  assert(eq_tc(lv4_lPpx,lv4_x));
  for i_0 in lv4keys do
    if (lv4_lPpx[i_0] ne 0) then
      mu_lpow:=lv4_x[i_0]/lv4_lPpx[i_0];
      break i_0;
    end if;
  end for;
  _<x>:=PolynomialRing(Parent(mu_lpow));
  mu:=RootsInSplittingField(x^l-mu_lpow)[1][1];
  for i in lv4keys do
    lv4_Ppx[i]*:=mu;
  end for;
  return lv4_P,lv4_Ppx;
end function;



//for x+e.
function calculate_nu_lpow(lv4tnp,lv4_e,lv4_x,lv4_xpe,l)
  lv4_lepx:=ThreePtLadder_plus(lv4tnp,l,lv4_e,lv4_x,lv4_xpe);
  for i_0 in lv4keys do
    if (lv4_lepx[i_0] ne 0) then
      nu_lpow:=lv4_x[i_0]/lv4_lepx[i_0];
      return nu_lpow;
    end if;
  end for;
end function;







//e in A[l],
//x in A.
//check if x+e,e are excellent for x.
function check_excellent_pt(lv4tnp,lv4tc_e,lv4tc_x,lv4tc_epx,l)
  return eq_Assoc(ThreePtLadder_plus(lv4tnp,l,lv4tc_e,lv4tc_x,lv4tc_epx),lv4tc_x);
end function;




//[LR22]Alg2.
//give excellent kernel K.
function excellent_Ker(lv4tnp,lv4tc_e1,lv4tc_e2,lv4tc_e12,l)
  time1:=Time();
  lv4tc_e1:=excellent_torsion(lv4tnp,lv4tc_e1,l);
  lv4tc_e2:=excellent_torsion(lv4tnp,lv4tc_e2,l);
  lv4tc_e12:=excellent_torsion(lv4tnp,lv4tc_e12,l);
  "Time excellent point.",Time(time1);
  
  time1:=Time();
  want:=linear_combination(lv4tnp,l-1,lv4tc_e1,lv4tc_e2,lv4tc_e12);
  "Time lin_com",Time(time1);

  return want;
end function;



//[LR22].Alg2.
//give excellent (x+K)^{~} for K^{~}.
function excellent_xpKer(lv4tnp,lv4tc_e1,lv4tc_e2,lv4tc_e12,lv4tc_x,lv4tc_e1px,lv4tc_e2px,l)

  time1:=Time();
  lv4tc_e12:=excellent_torsion(lv4tnp,lv4tc_e12,l);
  lv4tc_e1,lv4tc_e1px:=excellent_pt(lv4tnp,lv4tc_e1,lv4tc_x,lv4tc_e1px,l);
  lv4tc_e2,lv4tc_e2px:=excellent_pt(lv4tnp,lv4tc_e2,lv4tc_x,lv4tc_e2px,l);
  "Time excellent  x+e",Time(time1);


  time1:=Time();
  xpKer:=x_plus_lincom(lv4tnp,l-1,lv4tc_e1,lv4tc_e2,lv4tc_e12,lv4tc_x,lv4tc_e1px,lv4tc_e2px);
  "Time lincom x+K",Time(time1);

  return xpKer;
end function;




//give sum of square representation of l.
function sum_sq(l)
  assert(IsPrime(l));
  if l eq 2 then
    return [1,1];
  elif (l mod 4) eq 1 then
    for a in {1..l} do
      if IsSquare(l-a^2) then
        _,b:=IsSquare(l-a^2);
        assert(a^2+b^2 eq l);
        return [a,b];
      end if;
    end for;
  elif (l mod 4) eq 3 then
    for a in {1..l} do
      for b in {1..l} do
        for c in {1..l} do
          if IsSquare(l-a^2-b^2-c^2) then
            _,d:=IsSquare(l-a^2-b^2-c^2);
            return [a,b,c,d];
          end if;
        end for;
      end for;
    end for;
  end if;
end function;


        
/*

function level_up(lv4tnp,lv4tc_e1,lv4tc_e2,lv4tc_e12,lv4tc_x,lv4tc_e1px,lv4tc_e2px,l)
  assert(IsPrime(l));
  assert(l ne 2);
  sum_sq:=sum_sq(l);
  r:=#sum_sq;
  lv4lkeys:={[a,b]:a,b in {0..(4*l-1)}};
  lv4ltc:=AssociativeArray();
  l_torsion:=excellent_Ker(lv4tnp,lv4tc_e1,lv4tc_e2,lv4tc_e12,l)

  xpltor:=excellent_xpKer(lv4tnp,lv4tc_e1,lv4tc_e2,lv4tc_e12,lv4tc_x,lv4tc_e1px,lv4tc_e2px,l);

  for J in lv4lkeys do
    j:=[(J[1] mod 4),(J[2] mod 4)];
    q:=[(J[1] mod l),(J[2] mod l)];
    Q:=l_torsion[q[1],q[2]];


*/



function codomain_tnp(lv4tnp,lv4_e1,lv4_e2,lv4_e12,l)
  ex_K:=excellent_Ker(lv4tnp,lv4_e1,lv4_e2,lv4_e12,l);
  lv4tnp_cod:=AssociativeArray();
  for j in lv4keys do
    lv4tnp_cod[j]:=0;
    for key in Keys(ex_K) do
      lv4tnp_cod[j]+:=(ex_K[key][j])^l;
    end for;
  end for;
  return lv4tnp_cod,ex_K[[1,0]],ex_K[[0,1]],ex_K[[1,1]];
end function;



//not need field extension.
function codomain_tnp_2(lv4tnp,lv4_e1,lv4_e2,lv4_e12,l)
  time0:=Time();
  mu_1_lpow:=calculate_mu_lpow(lv4tnp,lv4_e1,l);
  mu_2_lpow:=calculate_mu_lpow(lv4tnp,lv4_e2,l);
  mu_12_lpow:=calculate_mu_lpow(lv4tnp,lv4_e12,l);
  "take lpow",Time(time0);
  time0:=Time();
  lincom_e1e2:=linear_combination(lv4tnp,l-1,lv4_e1,lv4_e2,lv4_e12);
  "calc lincom",Time(time0);
  lv4tnp_cod:=AssociativeArray();
  for j in lv4keys do
    lv4tnp_cod[j]:=0;
    for k1 in {0..(l-1)} do
      for k2 in {0..(l-1)} do
        lv4tnp_cod[j]+:=mu_1_lpow^(k1^2-k1*k2)*mu_2_lpow^(k2^2-k1*k2)*mu_12_lpow^(k1*k2)*(lincom_e1e2[[k1,k2]][j]^l);
      end for;
    end for;
  end for;
  return lv4tnp_cod,mu_1_lpow,mu_2_lpow,mu_12_lpow;
end function;








//[LR22]Alg4 for g=2,n=4.
function ll_isogeny(lv4tnp,lv4_e1,lv4_e2,lv4_e12,lv4_x,lv4_e1px,lv4_e2px,l)
  //first, construct excellent lift of K.

  time1:=Time();
  ex_xpK:=excellent_xpKer(lv4tnp,lv4_e1,lv4_e2,lv4_e12,lv4_x,lv4_e1px,lv4_e2px,l);
  "AllTime excellent x+ker",Time(time1);

  time1:=Time();
  lv4tc_fx:=AssociativeArray();
  for j in lv4keys do
    lv4tc_fx[j]:=0;
    for key in Keys(ex_xpK) do
      lv4tc_fx[j]+:=(ex_xpK[key][j])^l;
    end for;
  end for;
  "AllTime calc isogeny",Time(time1);

  return lv4tc_fx;
end function;



//speed up.
function ll_isogeny_2(lv4tnp,lv4_e1,lv4_e2,lv4_e12,lv4_x,lv4_e1px,lv4_e2px,l,mu_1_lpow,mu_2_lpow,mu_12_lpow)
  nu_1_lpow:=calculate_nu_lpow(lv4tnp,lv4_e1,lv4_x,lv4_e1px,l);
  nu_2_lpow:=calculate_nu_lpow(lv4tnp,lv4_e2,lv4_x,lv4_e2px,l);
  xpKer:=x_plus_lincom(lv4tnp,l-1,lv4_e1,lv4_e2,lv4_e12,lv4_x,lv4_e1px,lv4_e2px);
  lv4_fx:=AssociativeArray();
  for j in lv4keys do
    lv4_fx[j]:=0;
    for k1 in {0..(l-1)} do
      for k2 in {0..(l-1)} do
        excell_coff:=mu_1_lpow^(k1*(k1-k2-1))*mu_2_lpow^(k2*(k2-k1-1))*mu_12_lpow^(k1*k2)*nu_1_lpow^k1*nu_2_lpow^k2;
        lv4_fx[j]+:=excell_coff*(xpKer[[k1,k2]][j])^l;
      end for;
    end for;
  end for;
  return  lv4_fx;
end function;




