// FastThetasGenus2.m
//     by Hugo Labrande
//
// Distributed under a GPL license. See the paper
//         Hugo Labrande, Emmanuel Thomé, "Computing theta functions in quasi-linear time in genus 2 and above"
// for an explanation of the algorithm (or my PhD Thesis).
// If you use this code, please cite that paper (instead of just this file) in your paper. Thank you!


//=============== Not-so-naive algorithm ==================

function NaiveThetasGenus2(z,tau)
  // Fast-ish computation of theta : number of mults = 6 B + (8+4) sqrt(B) + a few, to compute all 4 thetas, and O(1) memory
  // Possible to get down to 2 or 3 multiplications if we allow O(B) memory

  // Assume 2 Im(tau_3) \leq Im(tau_1) \leq Im(tau_2)		(Minkowski-reduced)
  p := Precision(Parent(z[1][1]));
  RR := RealField(p);
  // Compute B
  B := RR!(2*p*Log(10)/Pi(RR)) + 3;		// Log(10) to get digits of precision, Log(2) to get bits
  // Get the precision right to counter rounding error
  //        in genus 1 this was p + 7log B, let's do p + 20log B, just in case, even if it's probably fine...
  CC := ComplexField(p+Ceiling(20*Log(B)));

  q1 := Exp(CC.1*Pi(CC)*tau[1][1]);	q1sq := q1^2;
  q2 := Exp(CC.1*Pi(CC)*tau[2][2]);	q2sq := q2^2;
  q3 := Exp(CC.1*Pi(CC)*tau[1][2]);	q3sq := q3^2;

  w1 := Exp(CC.1*Pi(CC)*z[1][1]);	w1sq := w1^2;
  w2 := Exp(CC.1*Pi(CC)*z[2][1]);	w2sq := w2^2;

  // 4 theta-constants with a=0
  thetaFunctions := [CC!1, CC!1, CC!1, CC!1];
  thetaConstants := [CC!1, CC!1, CC!1, CC!1];

  // n=0
  q1m2 := q1;
  q12mminus2 := q1sq;		// at the beginning of the loop, it contains 2(m+1)-2, because we use the relation alpha_m = f(alpha_(m-1), alpha_(m-2))
  r1 := w1sq + 1/w1sq; r1m := q1*r1; r1mminus1 := 2;
  whenToStop := Ceiling(Sqrt(B/Imaginary(tau[1][1]))) + 3;
  for m := 1 to whenToStop do
	// theta functions
	thetaFunctions[1] := thetaFunctions[1] +r1m; thetaFunctions[2] := thetaFunctions[2] +r1m;
	if (m mod 2 eq 1) then thetaFunctions[3] := thetaFunctions[3] -r1m; thetaFunctions[4] := thetaFunctions[4] -r1m; else thetaFunctions[3] := thetaFunctions[3] +r1m; thetaFunctions[4] := thetaFunctions[4] +r1m; end if;
	q12mminus2timesq1 := (q12mminus2*q1);	// mini-optimization
	bubu := r1m;
	r1m := r1*r1m*q12mminus2timesq1 - r1mminus1*(q12mminus2^2);
	r1mminus1 := bubu;
	// theta constants
	term := 2*q1m2;
	thetaConstants[1] := thetaConstants[1] +term; thetaConstants[2] := thetaConstants[2] +term;
	if (m mod 2 eq 1) then thetaConstants[3] := thetaConstants[3] -term; thetaConstants[4] := thetaConstants[4] -term; else thetaConstants[3] := thetaConstants[3] +term; thetaConstants[4] := thetaConstants[4] +term; end if;
	q1m2 := q1m2*q12mminus2timesq1;	// m^2 + 2m + 1 = (m+1)^2
	q12mminus2 := q12mminus2*q1sq;
  end for;

  // m,n >=1
  v1 := q3sq + 1/q3sq;
  vnminus1 := 2; vn := v1;
  q2to2nminus2 := 1;
  q2n2 := q2;

  q3to2n := q3sq;

  s1 := w2sq + 1/w2sq;

  // the smallest bottomleft square
  w1w2sq := w1sq*w2sq; w1invw2sq := w2sq/w1sq;
  q2s1 := q2*s1; q1r1 := q1*r1;  q1q2 := q1*q2;		// small optimisations
  betas := Matrix(2,2,[q1r1, q1q2*q3sq*(w1w2sq + 1/w1w2sq), CC!2, q2s1]);		// (1,0), (1,1), (0,0), (0,1)
  betaprimes := Matrix(2,2,[q1r1, (q1q2/q3sq)*(w1invw2sq + 1/w1invw2sq), CC!2, q2s1]);

  for n := 1 to Ceiling(Sqrt(B/Imaginary(tau[2][2])))+3 do
	// Compute the terms
	//printf "\n n = %o",n;
	if (n gt 3) then  whenToStop := B - (n-3)^2*Imaginary(tau[2][2]); else  whenToStop := B; end if;
		if (whenToStop le 0) then whenToStop:=0; end if;	whenToStop := Ceiling(Sqrt(whenToStop/Imaginary(tau[1][1]))) + 3;

	// thetaconstants
	q1to2mminus2 := q1sq;	// this squared gives q^(4m-4)
	alphamminus1 := 2*q2n2;
		// don't forget to add that to the sum (m=0)
		thetaConstants[1] := thetaConstants[1] +alphamminus1; thetaConstants[3] := thetaConstants[3] +alphamminus1;
		if (n mod 2 eq 1) then thetaConstants[2] := thetaConstants[2] -alphamminus1; thetaConstants[4] := thetaConstants[4] -alphamminus1;
		else thetaConstants[2] := thetaConstants[2] +alphamminus1; thetaConstants[4] := thetaConstants[4] +alphamminus1; end if;
	alpham := q1*q2n2*vn;

	// theta(z)
		term := betas[2][2];		// not betas+betaprimes (since m=0 we only add the term once);
		thetaFunctions[1] := thetaFunctions[1] + term;
		thetaFunctions[3] := thetaFunctions[3] + term;
		if (n mod 2 eq 1) then
			thetaFunctions[2] := thetaFunctions[2] - term; thetaFunctions[4] := thetaFunctions[4] - term;
		else 	thetaFunctions[2] := thetaFunctions[2] + term; thetaFunctions[4] := thetaFunctions[4] + term;	end if;

	alphazm := betas[1][2]; alphaprimezm := betaprimes[1][2];
	alphazmminus1 := betas[2][2]; alphaprimezmminus1 := betaprimes[2][2];

	for m:=1 to whenToStop do

           // thetaConstants
	   term := 2*alpham;
	   thetaConstants[1] := thetaConstants[1] + term;
           if (n mod 2 eq 1) then thetaConstants[2] := thetaConstants[2] -term; else thetaConstants[2] := thetaConstants[2] +term; end if;
           if (m mod 2 eq 1) then thetaConstants[3] := thetaConstants[3] -term; else thetaConstants[3] := thetaConstants[3] +term; end if;
           if ((m+n) mod 2 eq 1) then thetaConstants[4] := thetaConstants[4] -term; else thetaConstants[4] := thetaConstants[4] +term; end if;
	   // induction relation on alpham
	   bubu := alpham;
	   alpham := vn*(q1to2mminus2*q1)*alpham - (q1to2mminus2^2)*alphamminus1;
	   alphamminus1 := bubu;

	   // thetaFunctions
	   term := alphazm + alphaprimezm;
	   thetaFunctions[1] := thetaFunctions[1] + term;
	   if (n mod 2 eq 1) then thetaFunctions[2] := thetaFunctions[2] - term; else thetaFunctions[2] := thetaFunctions[2] + term; end if;
	   if (m mod 2 eq 1) then thetaFunctions[3] := thetaFunctions[3] - term; else thetaFunctions[3] := thetaFunctions[3] + term; end if;
	   if ((m+n) mod 2 eq 1) then thetaFunctions[4] := thetaFunctions[4] - term; else thetaFunctions[4] := thetaFunctions[4] + term; end if;
	   // induction relations
	   r1Xq1to2mminus2Xq1 := r1*q1to2mminus2*q1;
	   bubu := alphazm;
	   alphazm := alphazm*r1Xq1to2mminus2Xq1*q3to2n - (q1to2mminus2*q3to2n)^2*alphazmminus1;
	   alphazmminus1 := bubu;
	   bubu := alphaprimezm;
	   alphaprimezm := alphaprimezm*r1Xq1to2mminus2Xq1/q3to2n - (q1to2mminus2/q3to2n)^2*alphaprimezmminus1;
	   alphaprimezmminus1 := bubu;


	   // update the one that depend on m
	   q1to2mminus2 := q1to2mminus2*q1sq;
	end for;

	// finally update the ones that depend on n
	// thetaconstants
	bubu := vn;
	vn := vn*v1-vnminus1;
	vnminus1 := bubu;
	q2to2nminus2 := q2to2nminus2*q2sq;
	q2n2 := q2n2*q2to2nminus2*q2;
	// thetafunctions
	s1Xq2to2nminus2Xq2 := s1*q2to2nminus2*q2;
	bubu := [betas[1][2], betas[2][2]];
	betas[1][2] := betas[1][2]*s1Xq2to2nminus2Xq2*q3sq - betas[1][1]*(q2to2nminus2*q3sq)^2;
	betas[2][2] := betas[2][2]*s1Xq2to2nminus2Xq2 - q2to2nminus2^2*betas[2][1];
	betas[1][1] := bubu[1]; betas[2][1] := bubu[2];
	bubu := [betaprimes[1][2], betaprimes[2][2]];
	betaprimes[1][2] := betaprimes[1][2]*s1Xq2to2nminus2Xq2/q3sq - betaprimes[1][1]*(q2to2nminus2/q3sq)^2;
	betaprimes[2][2] := betaprimes[2][2]*s1Xq2to2nminus2Xq2 - q2to2nminus2^2*betaprimes[2][1];
	betaprimes[1][1] := bubu[1]; betaprimes[2][1] := bubu[2];
	q3to2n := q3to2n*q3sq;
	
  end for;

  Cfinal := ComplexField(p);
  thetaFunctions := [ Cfinal!thetaFunctions[i] : i in [1..4]];
  thetaConstants := [ Cfinal!thetaConstants[i] : i in [1..4]];
  return thetaFunctions, thetaConstants;
end function;

/*
prec := 1000;
C := ComplexField(prec);
val := C!0.866;
tau := Matrix(2,2, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, C!0.15+C.1*val/9, -C!0.4+C.1*val]);
z := Matrix(2,1,[C!val/5, -C!val/7]);


thetaFunctions, thetaConstants := NaiveThetasGenus2(z,tau);


ComplexField(10)!Abs(Theta(Matrix(4,1,[C!0,C!0, C!0,C!0]), ZeroMatrix(C,2,1), tau)-thetaConstants[1]);
ComplexField(10)!Abs(Theta(Matrix(4,1,[C!0,C!0, C!0,C!1/2]), ZeroMatrix(C,2,1), tau)-thetaConstants[2]);
ComplexField(10)!Abs(Theta(Matrix(4,1,[C!0,C!0, C!1/2,C!0]), ZeroMatrix(C,2,1), tau)-thetaConstants[3]);
ComplexField(10)!Abs(Theta(Matrix(4,1,[C!0,C!0, C!1/2,C!1/2]), ZeroMatrix(C,2,1), tau)-thetaConstants[4]);

ComplexField(10)!Abs(Theta(Matrix(4,1,[C!0,C!0, C!0,C!0]), z, tau)-thetaFunctions[1]);
ComplexField(10)!Abs(Theta(Matrix(4,1,[C!0,C!0, C!0,C!1/2]), z, tau)-thetaFunctions[2]);
ComplexField(10)!Abs(Theta(Matrix(4,1,[C!0,C!0, C!1/2,C!0]), z, tau)-thetaFunctions[3]);
ComplexField(10)!Abs(Theta(Matrix(4,1,[C!0,C!0, C!1/2,C!1/2]), z, tau)-thetaFunctions[4]);
*/




// =================== Our algorithm =====================





// ThetaGenus2(12,z,t) = Theta_12(z,t) & so on
function ThetaGenus2(n,z,t)
  Vec := [0/2,0/2,0/2,0/2];
  Vec[4] := Round(n mod 2)/2;
  n := Round((n - (n mod 2))/2);
  Vec[3] := Round(n mod 2)/2;
  n := Round((n - (n mod 2))/2);
  Vec[2] := Round(n mod 2)/2;
  n := Round((n - (n mod 2))/2);
  Vec[1] := Round(n mod 2)/2;
  
  return Theta( Matrix(4,1,Vec), z, t);
end function;



// SignTheta(n,z,t)*Sqrt(thetasquared/theta0squared) = theta/theta0
function SignTheta(n,z,t)
  Ctwenty := ComplexField(20);
  num := ThetaGenus2(n, Matrix(2,1,[Ctwenty!z[1][1], Ctwenty!z[2][1]]), Matrix(2,2, [Ctwenty!t[1][1], Ctwenty!t[1][2],Ctwenty!t[2][1],Ctwenty!t[2][2]]));
  den := ThetaGenus2(0, Matrix(2,1,[Ctwenty!z[1][1], Ctwenty!z[2][1]]), Matrix(2,2, [Ctwenty!t[1][1], Ctwenty!t[1][2],Ctwenty!t[2][1],Ctwenty!t[2][2]]));
  if (Real(num/den) lt 0) then return -1; else return 1; end if;
end function;




/*
prec := 60;
C := ComplexField(prec);
val := C!0.866;
tau := Matrix(2,2, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, C!0.15+C.1*val/9, -C!0.4+C.1*val]);
z := Matrix(2,1,[C!val/5, -C!val/7]);

//for i:=0 to 15 do
  //ComplexField(30)!ThetaGenus2(i,z,tau);
//end for;

a := Matrix(4,1, [ThetaGenus2(0,z,tau)^2, ThetaGenus2(1,z,tau)^2, ThetaGenus2(2,z,tau)^2, ThetaGenus2(3,z,tau)^2]);
b := Matrix(4,1, [ThetaGenus2(0,Matrix(2,1,[C!0,C!0]),tau)^2, ThetaGenus2(1,Matrix(2,1,[C!0,C!0]),tau)^2, ThetaGenus2(2,Matrix(2,1,[C!0,C!0]),tau)^2, ThetaGenus2(3,Matrix(2,1,[C!0,C!0]),tau)^2]);

sqrta := Matrix(4,1, [ThetaGenus2(0,z,tau), ThetaGenus2(1,z,tau), ThetaGenus2(2,z,tau), ThetaGenus2(3,z,tau)]);
sqrtb := Matrix(4,1, [ThetaGenus2(0,Matrix(2,1,[C!0,C!0]),tau), ThetaGenus2(1,Matrix(2,1,[C!0,C!0]),tau), ThetaGenus2(2,Matrix(2,1,[C!0,C!0]),tau), ThetaGenus2(3,Matrix(2,1,[C!0,C!0]),tau)]);
*/




function HadamardMatrix(fi, n)
  m := Matrix(2,2, [fi!1,fi!1,fi!1,fi!-1]);
  res := m;
  for i:=2 to n do
    res := TensorProduct(res,m);
  end for;
  return res;
end function;




// a = vector of 4 theta, b = vector of 4 theta constants
function F(a,b, z, t)
  n := #Rows(a);
  myC := ComplexField(22);
  // extract with the right sign
  rootA := Matrix(n, 1, [ Sqrt(a[i][1])*SignTheta(i-1, z, t) : i in [1..n]]);
  rootB := Matrix(n, 1, [ Sqrt(b[i][1])*SignTheta(i-1, ZeroMatrix(myC,2,1), t) : i in [1..n]]);

  // hadamard stuff for optimal computation
  hadam := HadamardMatrix(Parent(a[1][1]), 2);
  hadamardA := hadam*rootA;
  hadamardB := hadam*rootB;
  mults := Matrix(n,1,[hadamardA[i][1]*hadamardB[i][1] : i in [1..n]]);
  squares := Matrix(n,1,[hadamardB[i][1]^2 : i in [1..n]]);
  return hadam^(-1)*mults/4, hadam^(-1)*squares/4;
end function;

function Finfty(a,b,z,t)
  p := Precision(Parent(a[1][1]));
  myt := t;
  r := a; s := b;
  res := [[r,s]];
  while( Abs(s[1][1]-s[2][1]) gt 10^(-p+10) ) do
    r,s := F(r,s,z,myt);  Append(~res, [r,s]);
    myt := 2*myt;
  end while;
  r,s := F(r,s,z,myt);  Append(~res, [r,s]);
  return res;
end function;


function Pow2PowN(x, n)
  r := x;
  for i:=1 to n do r := r^2;  end for;
  return r;
end function;


function AGMPrime(a,b,z,t)
  R := Finfty(a,b,z,t);
  //S := [ [ (R[i][1]/R[i][3])^(2^(i-1)) * R[i][3],  (R[i][2]/R[i][4])^(2^(i-1)) * R[i][4] ] : i in [1..#R]];
  //lambda := S[#S][1];
  mu := R[#R][2][1][1];					// R[#R] = the last [r,s] ; R[#R][1] = 4,1 matrix
  qu := R[#R][1][1][1]/R[#R][2][1][1];
  //print "Calcul de lambda :";
  lambda := Pow2PowN(qu, #R-1) * R[#R][2][1][1];
  return lambda, mu;
end function;


function AllDuplication(a,b)
  // Given theta_{0,b}(z,t) and theta_{0,b}(0,t) compute theta_{a,b}(z,t)
  // Explicit formulas for genus 2 are found in Cosset

  // we still need to find the genus g way to do this (Makarov = sum ab, not sum (-1)ab ab)
  n := #a;
  ThetaProducts := Matrix(4,4, [a[1]*b[1],a[1]*b[2],a[1]*b[3],a[1]*b[4],
				a[2]*b[2],a[2]*b[1],a[2]*b[4],a[2]*b[3],
				a[3]*b[3],a[3]*b[4],a[3]*b[1],a[3]*b[2],
				a[4]*b[4],a[4]*b[3],a[4]*b[2],a[4]*b[1] ]);
  hadam := HadamardMatrix(Parent(a[1]), 2);
  ThetaProducts := hadam*ThetaProducts;
  ThetaProducts := ElementToSequence(ThetaProducts/4);

  // you need to add a minus sign in front of the odd theta-constants to get exactly cosset's formula
  ThetaProducts[5+1] := -ThetaProducts[5+1];
  ThetaProducts[7+1] := -ThetaProducts[7+1];
  ThetaProducts[10+1] := -ThetaProducts[10+1];
  ThetaProducts[11+1] := -ThetaProducts[11+1];
  ThetaProducts[13+1] := -ThetaProducts[13+1];
  ThetaProducts[14+1] := -ThetaProducts[14+1];

  return ThetaProducts;
end function;



function toInverse(a,b,z,t)
  // Given theta_i/theta_0 (z and 0), compute lambda1, lambda2, lambda3, mu1, mu2, mu3
  // ya ya, it's weird to have a function which goal is to compute z and tau but you give it z and tau (for the signs) in the args
  p := Precision(Parent(a[1][1]));
  n := #Rows(a);
  CC := ComplexField(p);

  // First step: compute 1/theta00(z)^2, 1/theta00(0)^2
  theta00z, theta000 := AGMPrime(a,b,z,t);
  theta00z := 1/theta00z; theta000 := 1/theta000;
  // then compute the other ones
  thetaZwithAequals0 := []; theta0withAequals0 := [];
  for i:=1 to 4 do Append(~thetaZwithAequals0, a[i][1]*theta00z); end for;
  for i:=1 to 4 do Append(~theta0withAequals0, b[i][1]*theta000); end for;

  //thetaZwithAequals0;
  //theta0withAequals0;

  // Then compute everything at 2tau (simpler conceptually and generalizable to genus g)
  rootA := [ Sqrt(thetaZwithAequals0[i])*SignTheta(i-1, z, t) : i in [1..n]];
  rootB := [ Sqrt(theta0withAequals0[i])*SignTheta(i-1, Matrix(2,1,[CC!0,CC!0]), t) : i in [1..n]];
  sixteenThetas := AllDuplication(rootA, rootB);
  sixteenThetaConstants := AllDuplication(rootB, rootB);
  
  // then give it to borchardt
  tt := 2*t;
  jm1 := Matrix(2,2, [-1-1/tt[1][1], -tt[1][2]/tt[1][1], -tt[1][2]/tt[1][1], tt[2][2]-tt[1][2]^2/tt[1][1] ]);
  jm2 := Matrix(2,2, [tt[1][1]-tt[2][1]^2/tt[2][2], -tt[1][2]/tt[2][2], -tt[1][2]/tt[2][2], -1-1/tt[2][2] ]);
  detdetdet := -Determinant(tt);
  jm3 := Matrix(2,2, [tt[1][1]/detdetdet, (tt[1][1]*tt[2][2]-tt[1][2]^2-tt[1][2])/detdetdet, (tt[1][1]*tt[2][2]-tt[1][2]^2-tt[1][2])/detdetdet, tt[2][2]/detdetdet  ]);

  zz1 := Matrix(2,1, [z[1][1]/tt[1][1],				z[1][1]*tt[1][2]/tt[1][1] - z[2][1]	]);
  zz2 := Matrix(2,1, [z[2][1]*tt[1][2]/tt[2][2] - z[1][1],	z[2][1]/tt[2][2]			]);
  zz3 := ZeroMatrix(CC, 2,1); // we don't care about z in the third case!

  // CAREFUL in Dupont's Phd he confuses M1 and M2 page 151 (tables switched) and page 197 
  lambda1, mu1 := AGMPrime( Matrix(4,1, [CC!1, sixteenThetas[9+1]/sixteenThetas[8+1], sixteenThetas[0+1]/sixteenThetas[8+1], sixteenThetas[1+1]/sixteenThetas[8+1]]),  Matrix(4,1, [CC!1, sixteenThetaConstants[9+1]/sixteenThetaConstants[8+1], sixteenThetaConstants[0+1]/sixteenThetaConstants[8+1], sixteenThetaConstants[1+1]/sixteenThetaConstants[8+1]]), zz1, jm1);

  lambda2, mu2 := AGMPrime( Matrix(4,1, [CC!1, sixteenThetas[0+1]/sixteenThetas[4+1], sixteenThetas[6+1]/sixteenThetas[4+1], sixteenThetas[2+1]/sixteenThetas[4+1]]),  Matrix(4,1, [CC!1, sixteenThetaConstants[0+1]/sixteenThetaConstants[4+1], sixteenThetaConstants[6+1]/sixteenThetaConstants[4+1], sixteenThetaConstants[2+1]/sixteenThetaConstants[4+1]]), zz2, jm2);

  lambda3, mu3 := AGMPrime( Matrix(4,1, [CC!1, sixteenThetas[8+1]/sixteenThetas[0+1], sixteenThetas[4+1]/sixteenThetas[0+1], sixteenThetas[12+1]/sixteenThetas[0+1]]), Matrix(4,1, [CC!1, sixteenThetaConstants[8+1]/sixteenThetaConstants[0+1], sixteenThetaConstants[4+1]/sixteenThetaConstants[0+1], sixteenThetaConstants[12+1]/sixteenThetaConstants[0+1]]), zz3, jm3);

  // then extract the stuff
  computedtau1 := CC.1/(mu1*sixteenThetaConstants[8+1]);
  computedtau2 := CC.1/(mu2*sixteenThetaConstants[4+1]);
  computedtau3 := 1/(mu3*sixteenThetaConstants[0+1]); //computedtau3 := Sqrt(computedtau3 + computedtau1*computedtau2);	// imaginary part positive (cause reduction)

  // because we used duplication formulas 
  computedtau1 := computedtau1/2; computedtau2 := computedtau2/2; computedtau3 := computedtau3/4;

  // Newton on lambda to avoid computing logs
  computedz1 := lambda1*(-CC.1*2*computedtau1*sixteenThetas[8+1]);
  computedz2 := lambda2*(-CC.1*2*computedtau2*sixteenThetas[4+1]);

  det2tau := 1/(-mu3*sixteenThetaConstants[0+1]);		// don't forget we're still working with 2tau at this point
  thirdformula := lambda3*det2tau*sixteenThetas[0+1];		// 2 x 2truc/4det(tau) = i pi


  return [1/computedz1, 1/computedz2, 1/thirdformula, computedtau1, computedtau2, computedtau3];
end function;


/*prec := 60;
C := ComplexField(prec);
val := C!0.866;
tau := Matrix(2,2, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, C!0.15+C.1*val/9, -C!0.4+C.1*val]);
z := Matrix(2,1,[C!val/5, -C!val/7]);
a := Matrix(4,1, [ThetaGenus2(0,z,tau)^2, ThetaGenus2(1,z,tau)^2, ThetaGenus2(2,z,tau)^2, ThetaGenus2(3,z,tau)^2]);
b := Matrix(4,1, [ThetaGenus2(0,Matrix(2,1,[C!0,C!0]),tau)^2, ThetaGenus2(1,Matrix(2,1,[C!0,C!0]),tau)^2, ThetaGenus2(2,Matrix(2,1,[C!0,C!0]),tau)^2, ThetaGenus2(3,Matrix(2,1,[C!0,C!0]),tau)^2]);

sqrta := Matrix(4,1, [ThetaGenus2(0,z,tau), ThetaGenus2(1,z,tau), ThetaGenus2(2,z,tau), ThetaGenus2(3,z,tau)]);
sqrtb := Matrix(4,1, [ThetaGenus2(0,Matrix(2,1,[C!0,C!0]),tau), ThetaGenus2(1,Matrix(2,1,[C!0,C!0]),tau), ThetaGenus2(2,Matrix(2,1,[C!0,C!0]),tau), ThetaGenus2(3,Matrix(2,1,[C!0,C!0]),tau)]);



quotientsZ := Matrix(4,1, [ a[i][1]/a[1][1] : i in [1..4]] );
quotients0 := Matrix(4,1, [ b[i][1]/b[1][1] : i in [1..4]] );

myres := toInverse(quotientsZ, quotients0, z, tau);

resultIwant := [ Exp(C.1*Pi(C)*z[1][1]^2/tau[1][1]), Exp(C.1*Pi(C)*z[2][1]^2/tau[2][2]), Exp(C.1*Pi(C)*((z[1][1]^2*tau[2][2] + z[2][1]^2*tau[1][1] - 2*z[1][1]*z[2][1]*tau[2][1])/Determinant(tau) - 1)), tau[1][1], tau[2][2], tau[1][2]^2-tau[1][1]*tau[2][2]];

[ RealField(10)!Abs(resultIwant[i]-myres[i]) : i in [1..6]];

tofind := (z[1][1]^2*tau[2][2] + z[2][1]^2*tau[1][1] - 2*z[1][1]*z[2][1]*tau[2][1])/Determinant(tau) - 1;
*/




function DiffFiniesOneStep(a,b, z,t, lambdaIwant, detIwant)
  // a is, as always, of size 4, but the first element is 1
  // Give approximations of the quotients (6 of them + 2 ones) with precision p and get an approximation with precision 2p
  prec := Precision(Parent(a[1][1]));
  CC := ComplexField(2*prec);

  newa := [ CC!a[i][1] : i in [1..#Rows(a)]];
  newb := [ CC!b[i][1] : i in [1..#Rows(b)]];


  myAofSize3 := [ newa[i]/newa[1] : i in [2..#newa]];
  myBofSize3 := [ newb[i]/newb[1] : i in [2..#newb]];

  epsilon := CC!(10^(-prec));
  // we're going to assume that z and t are still good enough to determine the correct sign
  toInverseDeBase := toInverse( Matrix(4,1,[CC!1,myAofSize3[1],myAofSize3[2],myAofSize3[3]]), Matrix(4,1,[CC!1,myBofSize3[1],myBofSize3[2],myBofSize3[3]]) ,z,t);

  perturb := [];
  Append(~perturb, toInverse( Matrix(4,1,[CC!1,myAofSize3[1]+epsilon,myAofSize3[2],myAofSize3[3]]), Matrix(4,1,[CC!1,myBofSize3[1],myBofSize3[2],myBofSize3[3]]) ,z,t));
  Append(~perturb, toInverse( Matrix(4,1,[CC!1,myAofSize3[1],myAofSize3[2]+epsilon,myAofSize3[3]]), Matrix(4,1,[CC!1,myBofSize3[1],myBofSize3[2],myBofSize3[3]]) ,z,t));
  Append(~perturb, toInverse( Matrix(4,1,[CC!1,myAofSize3[1],myAofSize3[2],myAofSize3[3]+epsilon]), Matrix(4,1,[CC!1,myBofSize3[1],myBofSize3[2],myBofSize3[3]]) ,z,t));
  Append(~perturb, toInverse( Matrix(4,1,[CC!1,myAofSize3[1],myAofSize3[2],myAofSize3[3]]), Matrix(4,1,[CC!1,myBofSize3[1]+epsilon,myBofSize3[2],myBofSize3[3]]) ,z,t));
  Append(~perturb, toInverse( Matrix(4,1,[CC!1,myAofSize3[1],myAofSize3[2],myAofSize3[3]]), Matrix(4,1,[CC!1,myBofSize3[1],myBofSize3[2]+epsilon,myBofSize3[3]]) ,z,t));
  Append(~perturb, toInverse( Matrix(4,1,[CC!1,myAofSize3[1],myAofSize3[2],myAofSize3[3]]), Matrix(4,1,[CC!1,myBofSize3[1],myBofSize3[2],myBofSize3[3]+epsilon]) ,z,t));

  // need to get an explicit form for the jacobian? bah, inverse(jacobian) is fine)
  jacobienne := ZeroMatrix(CC,6,6);
  for i:=1 to 6 do
    for j := 1 to 6 do
	jacobienne[i][j] := (perturb[j][i]-toInverseDeBase[i])/epsilon;
    end for;
  end for;
  //dispjac := jacobienne; dispjac := Matrix(6,6, [ ComplexField(3)!Abs(dispjac[i][j]) : i,j in [1..6]]); dispjac;

  changement := Matrix([[ toInverseDeBase[1]-lambdaIwant[1], toInverseDeBase[2]-lambdaIwant[2], toInverseDeBase[3]-lambdaIwant[3],
			  toInverseDeBase[4]-t[1][1], toInverseDeBase[5]-t[2][2], toInverseDeBase[6]-detIwant]]);

//  changement := Matrix([[ toInverseDeBase[1]-lambdaIwant[1], toInverseDeBase[2]-lambdaIwant[2], toInverseDeBase[3]-lambdaIwant[3],
//			  toInverseDeBase[4]-t[1][1], toInverseDeBase[5]-t[2][2], toInverseDeBase[6]-(t[1][2]^2-t[1][1]*t[2][2])]]);
  //RealField(10)!( Abs(changement[1][1]) + Abs(changement[1][2]) + Abs(changement[1][3]) + Abs(changement[1][4]) + Abs(changement[1][5]) + Abs(changement[1][6]) );
//  if ( RealField(10)!( Abs(changement[1][1]) + Abs(changement[1][2]) + Abs(changement[1][3]) + Abs(changement[1][4]) + Abs(changement[1][5]) + Abs(changement[1][6]) ) gt 10*Real(epsilon) ) then print "Something went wrong in the Newton"; end if;
  changement := changement*(Transpose(jacobienne)^(-1));
  //ComplexField(10)!( Abs(changement[1][1]) + Abs(changement[1][2]) + Abs(changement[1][3]) + Abs(changement[1][4]) + Abs(changement[1][5]) + Abs(changement[1][6]) );

  return Matrix(4,1,[1, myAofSize3[1]-changement[1][1], myAofSize3[2]-changement[1][2], myAofSize3[3]-changement[1][3]]), Matrix(4,1,[1, myBofSize3[1]-changement[1][4], myBofSize3[2]-changement[1][5], myBofSize3[3]-changement[1][6]]);
end function;


/*
prec := 600;
C := ComplexField(2*prec);
CC:= ComplexField(prec);
val := C!0.866;
tau := Matrix(2,2, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, C!0.15+C.1*val/9, -C!0.4+C.1*val]);
z := Matrix(2,1,[C!val/5, -C!val/7]);
val := CC!val;
tautau := Matrix(2,2, [CC!0.2+val*CC.1/2, CC!0.15+CC.1*val/9, CC!0.15+CC.1*val/9, -CC!0.4+CC.1*val]);
zz := Matrix(2,1,[CC!val/5, -CC!val/7]); 

//a := Matrix(4,1, [ThetaGenus2(0,z,tau)^2, ThetaGenus2(1,z,tau)^2, ThetaGenus2(2,z,tau)^2, ThetaGenus2(3,z,tau)^2]);
//b := Matrix(4,1, [ThetaGenus2(0,Matrix(2,1,[C!0,C!0]),tau)^2, ThetaGenus2(1,Matrix(2,1,[C!0,C!0]),tau)^2, ThetaGenus2(2,Matrix(2,1,[C!0,C!0]),tau)^2, ThetaGenus2(3,Matrix(2,1,[C!0,C!0]),tau)^2]);
a := Matrix(4,1, [ThetaGenus2(0,zz,tautau)^2, ThetaGenus2(1,zz,tautau)^2, ThetaGenus2(2,zz,tautau)^2, ThetaGenus2(3,zz,tautau)^2]);
b := Matrix(4,1, [ThetaGenus2(0,Matrix(2,1,[C!0,C!0]),tautau)^2, ThetaGenus2(1,Matrix(2,1,[C!0,C!0]),tautau)^2, ThetaGenus2(2,Matrix(2,1,[C!0,C!0]),tautau)^2, ThetaGenus2(3,Matrix(2,1,[C!0,C!0]),tautau)^2]);

lambdaIwant := [ Exp(C.1*Pi(C)*z[1][1]^2/tau[1][1]), Exp(C.1*Pi(C)*z[2][1]^2/tau[2][2]), Exp(C.1*Pi(C)*((z[1][1]^2*tau[2][2] + z[2][1]^2*tau[1][1] - 2*z[1][1]*z[2][1]*tau[2][1])/Determinant(tau) - 1))];

quotientsZ := Matrix(4,1, [ a[i][1]/a[1][1] : i in [1..4]] );
quotients0 := Matrix(4,1, [ b[i][1]/b[1][1] : i in [1..4]] );
//rrrrr := toInverse(quotientsZ, quotients0, z, tau);
// [ rrrrr[1]-z[1][1], rrrrr[2]-z[2][1], rrrrr[3]-tau[1][1], rrrrr[4]-tau[2][2], rrrrr[5]-tau[1][2], rrrrr[6]];

apresdiffA, apresdiffB := DiffFiniesOneStep(quotientsZ, quotients0, z, tau, lambdaIwant, tau[1][2]^2-tau[1][1]*tau[2][2]);



ComplexField(10)!apresdiffA[2][1]-ThetaGenus2(1,z,tau)^2/ThetaGenus2(0,z,tau)^2;
//C!quotientsZ[2][1]-ThetaGenus2(1,z,tau)^2/ThetaGenus2(0,z,tau)^2;

*/


LOW_PRECISION := 3000;

function CalculThetas(z, tau)

  // Returns [theta_{0,b}(z,tau)]_{b=0..3}, [theta_{0,b}(0,tau)]_{b=0...3}

  CC := Parent(z[1][1]);

  // determine which low precision (this is to avoid thresholding in the complexity)
  lowprec := Precision(CC);
  flag:=0;
  if (Precision(CC) lt LOW_PRECISION) then
	flag:=1;
  else
	while (lowprec gt LOW_PRECISION) do
	  lowprec := Round(lowprec/2) +10;		  // the "10" is to compensate precision losses, but it should be way too much
	end while;
  end if;

  // Low precision computation first
  Clow := ComplexField(lowprec);
  zerolow := ZeroMatrix(Clow,2,1);
  zlow := Matrix(2,1, [Clow!z[1][1], Clow!z[2][1]]);
  taulow := Matrix(2,2, [Clow!tau[1][1], Clow!tau[1][2], Clow!tau[2][1], Clow!tau[2][2]]);
  initA, initB := NaiveThetasGenus2(zlow,taulow);
  if (flag eq 1) then
	return initA, initB;
  else
	// we keep going
	initA := [ initA[i]^2 : i in [1..#initA]];   initB := [ initB[i]^2 : i in [1..#initB]];
  end if;


  // Computing the lamba, mu we're Newtonning on
  z1sq := z[1][1]^2;
  z2sq := z[2][1]^2;
  twoz1z2 := (z[1][1]+z[2][1])^2 - z1sq - z2sq;
  detIwant := tau[1][2]^2-tau[1][1]*tau[2][2];	// the factor that pops out for the third tau is -det
  IPI := CC.1*Pi(CC);
//  lambdaIwant := [ Exp(CC.1*Pi(CC)*z[1][1]^2/tau[1][1]), Exp(CC.1*Pi(CC)*z[2][1]^2/tau[2][2]), Exp(CC.1*Pi(CC)*((z[1][1]^2*tau[2][2] + z[2][1]^2*tau[1][1] - 2*z[1][1]*z[2][1]*tau[2][1])/Determinant(tau) - 1))];
  lambdaIwant := [ Exp(IPI*z1sq/tau[1][1]), Exp(IPI*z2sq/tau[2][2]), Exp(IPI*((z1sq*tau[2][2] + z2sq*tau[1][1] - twoz1z2*tau[2][1])/(-detIwant) - 1))];

  a := Matrix(4,1, [ initA[i]/initA[1] : i in [1..4]]);
  b := Matrix(4,1, [ initB[i]/initB[1] : i in [1..4]]);
  // Finite differences to refine the result
  p := lowprec;
  while (p lt Precision(z[1][1])) do
    p := 2*p;
    Cp := ComplexField(p);
    a, b := DiffFiniesOneStep(a,b, Matrix(2,1, [Cp!z[1][1], Cp!z[2][1]]), Matrix(2,2, [Cp!tau[1][1], Cp!tau[1][2], Cp!tau[2][1], Cp!tau[2][2]]), [Cp!lambdaIwant[1], Cp!lambdaIwant[2], Cp!lambdaIwant[3]], Cp!detIwant);
  end while;

  // We have the quotient of thetas, apply F to unstick them
  theta00z, theta000 := AGMPrime(a,b,z,tau);
  //theta00z, theta000;
  a := [ Sqrt(a[i][1]/theta00z)*SignTheta(2,z,tau) : i in [1..4]];
  b := [ Sqrt(b[i][1]/theta000)*SignTheta(2,ZeroMatrix(CC,2,1),tau) : i in [1..4]];

  a := [CC!a[1], CC!a[2], CC!a[3], CC!a[4]];
  b := [CC!b[1], CC!b[2], CC!b[3], CC!b[4]];
  return a,b;
end function;


function Magma8Thetas(z,tau)
  theC := ComplexField(Precision(z[1][1]));
  thezero := ZeroMatrix(theC,2,1);
  resA := [ThetaGenus2(0,z,tau), ThetaGenus2(1,z,tau), ThetaGenus2(2,z,tau), ThetaGenus2(3,z,tau)];
  resB := [ThetaGenus2(0,thezero,tau), ThetaGenus2(1,thezero,tau), ThetaGenus2(2,thezero,tau), ThetaGenus2(3,thezero,tau)];
  return resA, resB;
end function;


/*prec := 3200;
C := ComplexField(prec);
val := C!0.866;
tau := Matrix(2,2, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, C!0.15+C.1*val/9, -C!0.4+C.1*val]);
z := Matrix(2,1,[C!val/5, -C!val/7]);

time a,b := CalculThetas(z,tau);
time naiveA, naiveB := NaiveThetasGenus2(z,tau);
time expectedA, expectedB := Magma8Thetas(z,tau);
time notNeeded := ThetaGenus2(0,z,tau);

[ RealField(10)!Abs(a[i]-expectedA[i]) : i in [1..4]];
[ RealField(10)!Abs(b[i]-expectedB[i]) : i in [1..4]];
*/

for i:= 0 to 10 do
	prec := 1000*2^i;
        printf "\nPrec : %o\n",prec;
	C := ComplexField(prec);
	val := C!0.866;
	tau := Matrix(2,2, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, C!0.15+C.1*val/9, -C!0.4+C.1*val]);
	z := Matrix(2,1,[C!val/5, -C!val/7]);

	time a,b := CalculThetas(z,tau);
	time naiveA, naiveB := NaiveThetasGenus2(z,tau);
	//time expectedA, expectedB := Magma8Thetas(z,tau);
	time notNeeded := ThetaGenus2(0,z,tau);
end for;






















