//----------------------------------

procedure check_parameter(p,N_A,N_B)
  assert(IsPrime(p));
  "bit length of p:",Ilog2(p);
  assert(p mod 4 eq 3);
  assert(N_A gt N_B);
  assert(IsDivisibleBy(p+1,N_A) or IsDivisibleBy(p-1,N_A));
  assert(IsDivisibleBy(p+1,N_B) or IsDivisibleBy(p-1,N_B));
  assert(IsOdd(N_A));
  if (N_A-N_B)*N_B lt p then
    "we have to multiply.";
    if IsDivisibleBy(p+1,N_B) then
      pN_B:=p+1;
    else
      pN_B:=p-1;
    end if;
    fatoriztion_seq(pN_B);
  end if;
end procedure;


//----------------------------------
//5ビット安全.
p:=104959;//16bit
N_A:= 3^2 * 7^3 * 17;
N_B:=2^9 * 5;

//10ビット安全.
p:=202546499;//27bit
N_A:=7^2 * 11^2 * 19 * 29 * 31;
N_B:=2^2 * 3 * 5^3 * 13^2 * 17;

//20ビット安全.
p := 625750366823999; //49bit
N_A := 13 * 17^3 * 19^2 * 41 * 43^2;
N_B := 2^5 * 3^3 * 5^3 * 11 * 23 * 29 * 31 * 47;

//25ビット安全
p := 6510321409315018751;//62bit
N_A := 5^5 * 7^3 * 13 * 19 * 53 * 59 * 71 * 79;
N_B:= 2^11 * 3^2 * 17^3 * 29^2 * 37^2 * 41;

//----------------------------------






p:=991;//9bit
N_A:=3^2*5*11;
N_B:=2^5;

p:=911;//9bit
N_A:=5*7*13;
N_B:=2^4*3;

p:=859;//9bit
N_A:=3*11*13;
N_B:=2^2*5;

p:=11287;//13bit
N_A:=3^3*11*19;
N_B:=136;

p:=14951;//13bit
N_A:=5^2*13*23;
N_B:=168;


//小貫先生5ビット安全.
p:=104959;//16bit
N_A:= 3^2 * 7^3 * 17;
N_B:=2^9 * 5;


p:=1479871;//20bit
N_A:=3^6*5*7*29;
N_B:=1216;

//小貫先生10ビット安全.
p:=202546499;//27bit
N_A:=7^2 * 11^2 * 19 * 29 * 31;
N_B:=2^2 * 3 * 5^3 * 13^2 * 17;

//小貫先生20ビット安全.
p := 625750366823999;
N_A := 13 * 17^3 * 19^2 * 41 * 43^2;
N_B := 2^5 * 3^3 * 5^3 * 11 * 23 * 29 * 31 * 47;


//小貫先生25ビット安全.
p:=6510321409315018751;
N_A:=5^5 * 7^3 * 13 * 19 * 53 * 59 * 71 * 79;
N_B:=2^11 * 3^2 * 17^3 * 29^2 * 37^2 * 41;


//小貫先生30ビット安全.
p := 276154505650672190920223;
N_A:=13 * 23 * 37 * 43 * 47^2 * 59 * 61 * 71 * 73 * 97 * 101;
N_B:=2^5 * 3^6 * 7 * 11^4 * 19 * 29^3 * 67 * 79;




//小貫先生45ビット安全.
p := 69504748411397252246297776661471;
N_A := 5 * 7^2 * 19 * 41 * 59 * 101 * 137 * 157 * 227 * 233 * 317 * 347 * 349 * 457;
N_B := 2^5 * 3^7 * 11^4 * 17 * 29 * 97 * 149 * 163 * 257 * 373 * 401 * 439;





//B-SDIH論文example3改.
p:=0x76042798BBFB78AEBD02490BD2635DEC131ABFFFFFFFFFFFFFFFFFFFFFFFFFFF;//254bit
N_A:=3^(34)*11*17*19^2*53^2*97*107*109*131*137*197*199*227;
N_B:=2^110*5*7^2;


p:=0x78DAB06E306CA0903EF6085B501DF876D5BE579C27CE65FD5564603FBF88487F;//254bit
N_A:=3^2*5;
N_B:=2^5;


//SQISign論文.sectionC.
p:=2^32*5^21*7*11*163*1181*2389*5233*8353*10139*11939*22003*25391*41843*3726787*6548911-1;//256bit
N_A:=3^56;
N_B:=2^32*5^21*7*11;








//==================

//lに対する素数生成.
procedure compute_prime(min_p,max_p,max_l)
function smooth(p,max_l)
  facseq_p:=fatoriztion_seq(p+1);
  facseq_m:=fatoriztion_seq(p-1);
  Np:=1;
  Nm:=1;

  for i in {1..#facseq_p} do
    if facseq_p[i] le max_l then
      Np*:=facseq_p[i];
    end if;
  end for;

  for i in {1..#facseq_m} do
    if facseq_m[i] le max_l then
      Nm*:=facseq_m[i];
    end if;
  end for;

  N_A:=Max(Np,Nm);
  N_B:=1;

  if Np eq N_A then
    for i in {1..#facseq_m} do
      if N_B*facseq_m[i] le N_A then
        N_B*:=facseq_m[i];
      else
        break i;
      end if;
    end for;
  end if;
  
  if Nm eq N_A then
    for i in {1..#facseq_p} do
      if N_B*facseq_p[i] le N_A then
        N_B*:=facseq_p[i];
      else
        break i;
      end if;
    end for;
  end if;

  assert(N_A gt N_B);
  a:=N_A-N_B;
  if a*N_B le p then
    return false;
  end if;
  return true,N_A,N_B;
end function;

  for q in {min_p..max_p} do
    if IsPrime(q) and (q mod 4 eq 3) then
      if smooth(q,max_l) then
        _,N_A,N_B:=smooth(q,max_l);
        //if not(IsDivisibleBy(N_A,4)) then
          "prime",q,
          "N_A",fatoriztion_seq(N_A);
          N_B;
          break q;
          //fatoriztion_seq(N_B);
          "";
        //end if;
      end if;
    end if;
  end for;
end procedure;





//----------------------

function take_b(p,N_A,N_B)
  a:=N_A-N_B;
  assert(p mod 4 eq 3);
  assert(N_A gt N_B);
  if a*N_B gt p then
    times_FullRepInt_2(a*N_B,p,100);
    return "b",1;
  end if;
  if IsDivisibleBy(p+1,N_B) then
    ppm:=p+1;
  else
    ppm:=p-1;
  end if;
  assert(IsDivisibleBy(ppm,N_B));
  ppmdivNB:=(ppm div N_B);
  //Divisors((ppm div N_B));
  for b in Divisors((ppm div N_B)) do
    if (a*b*N_B gt p) then
      if GCD(b,N_A) eq 1 and GCD(b,a) eq 1 then
        times_FullRepInt_2(a*b*N_B,p,100);
        return "b",b;
      end if;
    end if;
  end for;
  return "Nothing.";
end function;



//----------------------
a:=N_A-N_B;
if IsDivisibleBy(p+1,N_B) then
  pN_B:=p+1;
else
  pN_B:=p-1;
end if;
b:=pN_B div N_B;
assert(IsDivisibleBy(pN_B,N_B));
assert(b*a*N_B gt p);
assert(b*N_B ge pN_B);
assert(GCD(b,N_A) eq 1);
assert(GCD(b,a) eq 1);
//--------------------



//------------------


function estimate(N_A)
  seq:=fatoriztion_seq(N_A);
  r:=[];
  for i in {1..#seq} do
    if seq[i] mod 4 eq 1 then
      r[i]:=2;
    elif seq[i] mod 4 eq 3 then
      r[i]:=4;
    else
      assert(false);
    end if;
  end for;
  sum_time:=0;
  for i in {1..#seq} do
    sum_time+:=(1500*(seq[i]^2+seq[i]^(5/2))+2*(seq[i]^r[i]));
  end for;
  return sum_time/(10^3);
end function;










