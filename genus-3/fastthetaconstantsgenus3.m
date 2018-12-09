// FastThetaConstantsGenus3.m
//     by Hugo Labrande
// (Version: 11/01/2018)
//
// Distributed under a GPL license.
// This is basically the same method as in Régis Dupont's PhD thesis, but with improvements/tricks given in
//         Hugo Labrande, Emmanuel Thomé, "Computing theta functions in quasi-linear time in genus 2 and above".
// Thanks to Christophe Ritzenthaler for testing the code and providing feedback.
// Thanks to Marco Streng for bug reports, and for his lemma on Minkowski-reduced matrices (used in the naive algo)
// Thanks to Enea Milio for comments and bug reports

// Documentation :
/*
  The numbering used all through the code is as follows:
    theta_i = theta [a, b] with i the integer whose binary expansion is the same as "ab"; this is Dupont's and Cosset's numbering
    ex: theta[ 0 1/2 0, 0 1/2 1/2] = theta_19
  
  Apologies in advance for the French words popping up in variable/function names.

  ** Use this code on period matrices that are Minkowski-reduced **
        This is needed for the first approximation (via the naive algorithm) to be accurate.
        
  ======== Functions

  CalculThetas is the main function; give it a period matrix tau, and it computes theta_i(0,\tau) for i=0...7.
    It uses the naive algorithm for precisions smaller than 450, then uses Newton's method to refine the approximation.
  NaiveThetaConstantsGenus3 is the naive algorithm; takes tau, and an optional value for some tests that are currently commented out (so put any value)
  ThetaGenus3 calls Magma's Theta function; give it the integer i and the period matrix tau

  ========= Computing all the thetas

  CalculThetas and NaiveThetaConstantsGenus3 only give the fundamental theta constants. To get all 64 theta constants, use the AllDuplication function, which applies the duplication formulas to give theta[a, b](0,2\tau)^2. Something like:

  // computation of theta[0,b](tau/2)
  fundamentalThetas := CalculThetas(tau/2);
  // computation of theta[a,b](tau)^2
  thetas := AllDuplication(fundamentalThetas);
  // square root extraction
  Clow := ComplexField(20); taulow := [ Clow!tau[i][j] : i,j in [1..3]];
  for i:=1 to 64 do
    thetas[i] := Sqrt(thetas[i]); if ( Real(ThetaGenus3(i,tau)) lt 0) then thetas[i] := -thetas[i]; end if;
  end for;
  
  ========= Refining the results

  If the results you get are not precise enough for your applications, you can refine them using Newton's method; this is the function DiffFiniesOneStep (as in 'one step of finite differences-based Newton').
  The arguments it takes are: the vector [1, theta_1^2/theta_0^2, theta_2^2/theta_0^2, ..., theta_7^2/theta_0^2], tau, and the output of ComputeTableOfSigns(tau) (precomputation of the signs of quotients of thetas, needed to extract correct square roots).

  // if b = [theta_0, theta_1, theta_2, ..., theta_7]
  quotients0 := Matrix(8,1, [ (b[i]/b[1])^2 : i in [1..8]] );
  apresdiffB := DiffFiniesOneStep(quotients0, tau, ComputeTableOfSigns(tau));

  If you have something accurate with P digits, the output will be squares of quotients with roughly 2P-5 digits of accuracy.
  
  ========= Timings

  Keep in mind that both algorithms compared here output 8 fundamental theta-constants.
  On my machine at work :
    precision 1000 (base-10 digits) : 5.7s vs 16.4 for naive algorithm (magma: 13.9 for a single theta-constant)
    precision 2000: 14.9s vs 117 for naive algorithm (magma: 117 for a single theta-constant)
    precision 4000: 46.9s vs 817 for naive algorithm (magma: 1145 just for one theta constant)
    precision 8000: 148s (others DNF)

  Happy genus 3 computing!
*/




function Terme(tau, i,j,k)
  CC := ComplexField(Precision(tau[1][1]));
  return Exp(CC.1*Pi(CC)*(i^2*tau[1][1] + j^2*tau[2][2] + k^2*tau[3][3] + 2*i*j*tau[1][2] + 2*i*k*tau[1][3] + 2*j*k*tau[2][3] ));
end function;


function SuperNaiveThetaConstantsGenus3(tau, bou)
  // The most naive function possible
  sum := 0;
  for i:=-bou to bou do
    for j:=-bou to bou do
      for k:=-bou to bou do
        sum := sum+Terme(tau,i,j,k);
      end for;
    end for;
  end for;
  return sum;
end function;




//=============== Not-so-naive algorithm ==================


function NaiveThetaConstantsGenus3(tau, check)
  // Fast-ish computation of theta : number of mults = 6 B + (8+4) sqrt(B) + a few, to compute all 4 thetas, and O(1) memory

  // Assume Minkowski-reduced
  assert ( Abs(tau[1][1]) gt 1 and Abs(tau[2][2]) gt 1 and Abs(tau[3][3]) gt 1
    and Abs(Imaginary(tau[1][1])) lt Abs(Imaginary(tau[2][2])) and Abs(Imaginary(tau[2][2])) lt Abs(Imaginary(tau[3][3]))
    and 2*Imaginary(tau[1][2]) lt Imaginary(tau[1][1]) and 2*Imaginary(tau[2][3]) lt Imaginary(tau[2][2]) or 2*Abs(Imaginary(tau[1][3])) lt Imaginary(tau[1][1]));

  
  p := Precision(Parent(tau[1][1]));
  RR := RealField(p);

  // Compute B
  //    we're summing on [-B,B]^3 (cube)
  //           TODO we could speed up the algorithm by finding an analysis for summing over the ellipsoid
  //    according to our analysis, theta-S_B < 24 e^{-pi c B^2}{(1-e^{-pi c})^3}
  c := Min( Min(Imaginary(tau[1][1]-tau[1][2]-Abs(tau[1][3])), Imaginary(tau[2][2]-tau[1][2]-tau[2][3])), Imaginary(tau[3][3]-tau[2][3]-Abs(tau[1][3])) );
  c := Max(c, Imaginary(tau[1][1])/100);
  B := RR!(Real(   (p+2-3*Log(1-Exp(-Pi(ComplexField(10))*c))/Log(10)) / (c*Pi(ComplexField(10))/Log(10)) ));

  //     we don't have a proof here, so this is experimental
  //     but we know we have bounds like e^(-Pi(C)*B^2)), so we look for c sqrt(P/pi)
  //     the Log(10) is there to get digits of precision (change to Log(2) for bits)
  //B := RR!((5/2)*p*Log(10)/Pi(RR)) + 3;		// Log(10) to get digits of precision, Log(2) to get bits
  B := Round(Sqrt(B))+1;
  //B:=5;
  //printf "B = %o\n",B;
  // Get the precision right to counter rounding error
  //        in genus 1 this was p + 7log B, let's do p + 20log B, just in case, even if it's probably fine...
  CC := ComplexField(p+Ceiling(20*Log(B)));
  ipi := CC.1*Pi(CC);

  q1 := Exp(ipi*tau[1][1]);	q1sq := q1^2;
  q2 := Exp(ipi*tau[2][2]);	q2sq := q2^2;
  q3 := Exp(ipi*tau[3][3]);	q3sq := q3^2;
  q4 := Exp(ipi*tau[1][2]);	q4sq := q4^2;  q4Invsq := 1/q4sq;
  q5 := Exp(ipi*tau[2][3]);	q5sq := q5^2;  q5Invsq := 1/q5sq;
  q6 := Exp(ipi*tau[1][3]);	q6sq := q6^2;  q6Invsq := 1/q6sq;

  // 8 theta-constants with a=0
  thetaConstants := [CC!1, CC!1, CC!1, CC!1, CC!1, CC!1, CC!1, CC!1];

  termI := CC!1;  // termI[1] : j=0, k=0, i>0 + i<0
  termJ := [CC!1,CC!1];  // termJ[1] : k=0, i>0, j>0 + i<0 j<0 ; termJ[2] : k=0, i>0, j<0 + i<0 j>0
  term := [CC!1,CC!1,CC!1,CC!1]; // term[1] : +++ & ---, term[2] : ++- & --+, term[3] : +-+ & -+-, term[4] : +-- & -++

  q12ip1 := q1;
  q42i := CC!1;		q4Inv2i := CC!1;
  q62i := CC!1;		q6Inv2i := CC!1;
  for i := 1 to B do
	//printf "i = %o \n",i;
    q22jp1 := q2;
    q52j := CC!1;	q5Inv2j := CC!1;
    termJ := [termI, termI];
    for j := 1 to B do			// TODO: this bound could be improved (ellipsoid)
      q32kp1 := q3;
      term := [termJ[1], termJ[1], termJ[2], termJ[2]];
      for k := 1 to B do		// TODO: this bound could be improved (ellipsoid)
        mytemp1 := q32kp1*q62i; mytemp2 := q32kp1*q6Inv2i;
        term[1] := term[1]*mytemp1*q52j;
        term[2] := term[2]*mytemp2*q5Inv2j;
        term[3] := term[3]*mytemp1*q5Inv2j;
        term[4] := term[4]*mytemp2*q52j;
        q32kp1 := q32kp1*q3sq;

	// Adding (i-1,j-1,k) (i-1,j-1,-k) (i-1,-(j-1),k) (i-1,-(j-1),-k)
	if (j eq 1) then
		toadd := term[1]+term[2];
	       	thetaConstants[1] := thetaConstants[1] + toadd;
		thetaConstants[3] := thetaConstants[3] + toadd;
		if (k mod 2 eq 1) then thetaConstants[2] := thetaConstants[2] - toadd; thetaConstants[4] := thetaConstants[4] - toadd; else thetaConstants[2] := thetaConstants[2] + toadd; thetaConstants[4] := thetaConstants[4] + toadd; end if;
		if ( (i-1) mod 2 eq 1) then thetaConstants[5] := thetaConstants[5] - toadd; thetaConstants[7] := thetaConstants[7] - toadd; else   thetaConstants[5] := thetaConstants[5] + toadd; thetaConstants[7] := thetaConstants[7] + toadd; end if;
		if ( (i-1+k) mod 2 eq 1) then thetaConstants[6] := thetaConstants[6] - toadd; thetaConstants[8] := thetaConstants[8] - toadd; else   thetaConstants[6] := thetaConstants[6] + toadd; thetaConstants[8] := thetaConstants[8] + toadd; end if;
	else
		toadd := term[1]+term[2]+term[3]+term[4];
       		thetaConstants[1] := thetaConstants[1] + toadd;
		if (k mod 2 eq 1) then thetaConstants[2] := thetaConstants[2]-toadd; else thetaConstants[2] := thetaConstants[2]+toadd; end if;
		if ( (j-1) mod 2 eq 1) then thetaConstants[3] := thetaConstants[3]-toadd; else thetaConstants[3] := thetaConstants[3]+toadd; end if;
		if ( (j-1+k) mod 2 eq 1) then thetaConstants[4] := thetaConstants[4]-toadd; else thetaConstants[4] := thetaConstants[4]+toadd; end if;
		if ( (i-1) mod 2 eq 1) then thetaConstants[5] := thetaConstants[5]-toadd; else thetaConstants[5] := thetaConstants[5]+toadd; end if;
		if ( (i-1+k) mod 2 eq 1) then thetaConstants[6] := thetaConstants[6]-toadd; else thetaConstants[6] := thetaConstants[6]+toadd; end if;
		if ( (i+j) mod 2 eq 1) then thetaConstants[7] := thetaConstants[7]-toadd; else thetaConstants[7] := thetaConstants[7]+toadd; end if;
		if ( (i+j+k) mod 2 eq 1) then thetaConstants[8] := thetaConstants[8]-toadd; else thetaConstants[8] := thetaConstants[8]+toadd; end if;
	end if;
	//printf "Diff (i-1,j-1,k) = (%o,%o,%o): %o\n",i-1,j-1,k,ComplexField(10)!Abs(thetaConstants[3]-check);

      end for;
      termJ[1] := termJ[1]*q22jp1*q42i;	// now it contains the term (i-1,j,0)
      termJ[2] := termJ[2]*q22jp1*q4Inv2i;     // now it contains the term (i-1,-j,0)

      // Adding (i-1,j,0) (i-1,-j,0)
	toadd := termJ[1]+termJ[2];
      thetaConstants[1] := thetaConstants[1] + toadd; thetaConstants[2] := thetaConstants[2] + toadd;
	if ( j mod 2 eq 1) then thetaConstants[3] := thetaConstants[3] - toadd; thetaConstants[4] := thetaConstants[4] - toadd; else thetaConstants[3] := thetaConstants[3] + toadd; thetaConstants[4] := thetaConstants[4] + toadd; end if;
	if ( (i-1) mod 2 eq 1) then thetaConstants[5] := thetaConstants[5] - toadd; thetaConstants[6] := thetaConstants[6] - toadd; else thetaConstants[5] := thetaConstants[5] + toadd; thetaConstants[6] := thetaConstants[6] + toadd; end if;
	if ( (i-1+j) mod 2 eq 1) then thetaConstants[7] := thetaConstants[7] - toadd; thetaConstants[8] := thetaConstants[8] - toadd; else thetaConstants[7] := thetaConstants[7] + toadd; thetaConstants[8] := thetaConstants[8] + toadd; end if;

	//printf "Diff (i-1,j,k) = (%o,%o,ZERO): %o\n",i-1,j,ComplexField(10)!Abs(thetaConstants[3]-check);

      q22jp1 := q22jp1*q2sq;
      q52j := q52j*q5sq;	q5Inv2j := q5Inv2j*q5Invsq;
    end for;

    termI := termI*q12ip1;
    q12ip1 := q12ip1*q1sq;
    q42i := q42i*q4sq;		q4Inv2i := q4Inv2i*q4Invsq;
    q62i := q62i*q6sq;		q6Inv2i := q6Inv2i*q6Invsq;

    if (i eq 1) then
      // i-1>1 from then on so -1 and 1 are different => term symmetry gives a factor 2
      termI := termI*2;
    end if;
    // Adding the term (i,0,0)
    thetaConstants[1] := thetaConstants[1] + termI;
    thetaConstants[2] := thetaConstants[2] + termI;
    thetaConstants[3] := thetaConstants[3] + termI;
    thetaConstants[4] := thetaConstants[4] + termI;
    if (i mod 2 eq 1) then thetaConstants[5] := thetaConstants[5] - termI; thetaConstants[6] := thetaConstants[6] - termI; thetaConstants[7] := thetaConstants[7] - termI; thetaConstants[8] := thetaConstants[8] - termI; else thetaConstants[5] := thetaConstants[5] + termI; thetaConstants[6] := thetaConstants[6] + termI; thetaConstants[7] := thetaConstants[7] + termI; thetaConstants[8] := thetaConstants[8] + termI; end if;
	//printf "Diff (i,j,k) = (%o,ZERO,ZERO): %o\n",i,ComplexField(10)!Abs(thetaConstants[3]-check);

  end for;


  Cfinal := ComplexField(p);
  thetaConstants := [ Cfinal!thetaConstants[i] : i in [1..8]];
  return thetaConstants;
end function;

/*
function NaiveThetaConstantsGenus3(tau)
  return NaiveThetaConstantsGenus3(tau, 1);
end function;
*/

// ThetaGenus3(12,t) = Theta_12(0,t) & so on
function ThetaGenus3(n,t)
  Vec := [0/2,0/2,0/2,0/2,0/2,0/2];
  Vec[6] := Round(n mod 2)/2;
  n := Round((n - (n mod 2))/2);
  Vec[5] := Round(n mod 2)/2;
  n := Round((n - (n mod 2))/2);
  Vec[4] := Round(n mod 2)/2;
  n := Round((n - (n mod 2))/2);
  Vec[3] := Round(n mod 2)/2;
  n := Round((n - (n mod 2))/2);
  Vec[2] := Round(n mod 2)/2;
  n := Round((n - (n mod 2))/2);
  Vec[1] := Round(n mod 2)/2; 
C := Parent(t[1][1]);
  return Theta( Matrix(6,1,Vec), ZeroMatrix(C,3,1), t);
end function;

/*
prec := 500;
C := ComplexField(prec);
val := C!0.866;
tau := Matrix(3,3, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, -C!0.3+C.1*val/5, C!0.15+C.1*val/9, -C!0.4+C.1*val, C.1*val/3, -C!0.3+C.1*val/5, C.1*val/3, C!0.05+C.1*(3*val/2)]);

time thetamagma := ThetaGenus3(0, tau);
//time thetasupernaive := SuperNaiveThetaConstantsGenus3(tau, 50);
//time thetaConstants := NaiveThetaConstantsGenus3(tau, thetasupernaive);
time thetaConstants := NaiveThetaConstantsGenus3(tau, ThetaGenus3(2,tau));
//time thetaConstants := NaiveThetaConstantsGenus3(tau, C!1);


ComplexField(10)!Abs(thetamagma-thetaConstants[1]);
//ComplexField(10)!Abs(thetasupernaive-thetaConstants[1]);



ComplexField(10)!Abs(ThetaGenus3(1, tau)-thetaConstants[2]);
ComplexField(10)!Abs(ThetaGenus3(2, tau)-thetaConstants[3]);
ComplexField(10)!Abs(ThetaGenus3(3, tau)-thetaConstants[4]);

ComplexField(10)!Abs(ThetaGenus3(4, tau)-thetaConstants[5]);
ComplexField(10)!Abs(ThetaGenus3(5, tau)-thetaConstants[6]);
ComplexField(10)!Abs(ThetaGenus3(6, tau)-thetaConstants[7]);
ComplexField(10)!Abs(ThetaGenus3(7, tau)-thetaConstants[8]);

*/




// =================== Our algorithm =====================




// Old version: recompute every time
//     this is wasteful and accounts for a big part of the computation
//     (ex with prec 400 : 35s with, 1.6s with the new version)
// SignTheta(n,t)*Sqrt(thetasquared/theta0squared) = theta/theta0

function SignTheta(n,t)
  Clow := ComplexField(10);
  taulow := Matrix(3,3, [ Clow!t[i][j] : i,j in [1..3]]);
  num := ThetaGenus3(n, taulow);
  den := ThetaGenus3(0, taulow);
  if (Real(num/den) lt 0) then return -1; else return 1; end if;
end function;


// New version: precompute the signs and give them in the functions
function ComputeOneTable(tau)
  Clow := ComplexField(10);
  taulow := Matrix(3,3, [ Clow!tau[i][j] : i,j in [1..3]]);
  mat := ElementToSequence(ZeroMatrix(Integers(),64,1));
  den := ThetaGenus3(0, taulow);
  for i:=0 to 7 do
    num := ThetaGenus3(i, taulow);
    if (Real(num/den) lt 0) then mat[i+1] := -1; else mat[i+1] := 1; end if;
  end for;
  return mat;
end function;

function ComputeTableOfSigns(tau)
  Clow := ComplexField(10);
  taulow := Matrix(3,3, [ Clow!tau[i][j] : i,j in [1..3]]);
  tt := 2*taulow;

  i3 := IdentityMatrix(Clow,3);
  delta11 := Matrix(3,3, [Clow!1,Clow!0,Clow!0,Clow!0,Clow!0,Clow!0,Clow!0,Clow!0,Clow!0]);
  delta22 := Matrix(3,3, [Clow!0,Clow!0,Clow!0,Clow!0,Clow!1,Clow!0,Clow!0,Clow!0,Clow!0]);
  delta33 := Matrix(3,3, [Clow!0,Clow!0,Clow!0,Clow!0,Clow!0,Clow!0,Clow!0,Clow!0,Clow!1]);
  delta12 := Matrix(3,3, [Clow!0,Clow!1,Clow!0,Clow!1,Clow!0,Clow!0,Clow!0,Clow!0,Clow!0]);
  delta13 := Matrix(3,3, [Clow!0,Clow!0,Clow!1,Clow!0,Clow!0,Clow!0,Clow!1,Clow!0,Clow!0]);
  delta23 := Matrix(3,3, [Clow!0,Clow!0,Clow!0,Clow!0,Clow!0,Clow!1,Clow!0,Clow!1,Clow!0]);
  jm1 := (-tt-delta11)*(delta11*tt+delta11-i3)^(-1);
  jm2 := (-tt-delta22)*(delta22*tt+delta22-i3)^(-1);
  jm3 := (-tt-delta33)*(delta33*tt+delta33-i3)^(-1);
  jm12 := (-tt-delta12)*(delta12*tt+delta11+delta22-i3)^(-1);
  jm13 := (-tt-delta13)*(delta13*tt+delta11+delta33-i3)^(-1);
  jm23 := (-tt-delta23)*(delta23*tt+delta33+delta22-i3)^(-1);
  ttinv := -(tau)^(-1);

  return [ ComputeOneTable(taulow), ComputeOneTable(jm1), ComputeOneTable(jm2), ComputeOneTable(jm3), ComputeOneTable(jm12), ComputeOneTable(jm13), ComputeOneTable(jm23), ComputeOneTable(ttinv)];
end function;

// We need this function because we can't precompute the tables for 2^k tau
//   Lemma (not proven... but it must be true...): if Re(theta_i(tau)/theta_0(tau)) > 0, Re(theta_i(2tau)/theta_0(2tau)) > 0 too
function GetCorrectSign(table, i, tau)
  if (table[i] eq 1) then return 1;
  else return SignTheta(i,tau);
  end if;
end function;


/*
prec := 60;
C := ComplexField(prec);
val := C!0.866;
tau := Matrix(3,3, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, C!0.05+C.1*val/10, C!0.15+C.1*val/9, -C!0.4+C.1*val, C!0.2+C.1*val/20, C!0.05+C.1*val/10, C!0.2+C.1*val/20, C!0.1+C.1*val*C!1.5]);

//for i:=0 to 15 do
  //ComplexField(30)!ThetaGenus2(i,z,tau);
//end for;

b := Matrix(4,1, [ThetaGenus3(0,tau)^2, ThetaGenus3(1,tau)^2, ThetaGenus3(2,tau)^2, ThetaGenus3(3,tau)^2]);

sqrtb := Matrix(4,1, [ThetaGenus3(0,tau), ThetaGenus3(1,tau), ThetaGenus3(2,tau), ThetaGenus3(3,tau)]);

table := ComputeOneTable(tau);

 Abs(Sqrt(b[1][1])*table[1] - sqrtb[1][1]);
 Abs(Sqrt(b[4][1])*table[4] - sqrtb[4][1]);
*/




function HadamardMatrix(fi, n)
  m := Matrix(2,2, [fi!1,fi!1,fi!1,fi!-1]);
  res := m;
  for i:=2 to n do
    res := TensorProduct(res,m);
  end for;
  return res;
end function;




// b = vector of 8 theta constants
function F(b, t, shortTableOfSigns)
  n := #Rows(b);
  // extract with the right sign
  rootB := Matrix(n, 1, [ Sqrt(b[i][1])*GetCorrectSign(shortTableOfSigns,i,t) : i in [1..n]]);
  // hadamard stuff for optimal computation
  hadam := HadamardMatrix(Parent(b[1][1]), 3);
  hadamardB := hadam*rootB;
  squares := Matrix(n,1,[hadamardB[i][1]^2 : i in [1..n]]);
  return hadam^(-1)*squares/n;
end function;

/*
prec := 400;
C := ComplexField(prec);
val := C!0.866;
tau := Matrix(3,3, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, C!0.05+C.1*val/10, C!0.15+C.1*val/9, -C!0.4+C.1*val, C!0.2+C.1*val/20, C!0.05+C.1*val/10, C!0.2+C.1*val/20, C!0.1+C.1*val*C!1.5]);
thconstants := Matrix(8, 1, [ ThetaGenus3(i,tau)^2 : i in [0..7]]);
time thconstants2tau := F(thconstants, tau, ComputeOneTable(tau));
ComplexField(10)!Abs(thconstants2tau[1][1]-ThetaGenus3(0,2*tau)^2);
*/

function Finfty(b,t, shortTableOfSigns)
  p := Precision(Parent(b[1][1]));
  n := #Rows(b);
  myt := t;
  s := b;
  res := [[s]];
  while( Abs(s[1][1]-s[n][1]) gt 10^(-p+10) ) do
    s := F(s,myt, shortTableOfSigns);  Append(~res, [s]);	// we use the table of signs as a first approximations of the signs at 2^k tau, using GetCorrectSign to correct this ; since we hope they're always all positive, it's not a bit hit in performance to have to recompute when it's negative
    myt := 2*myt;
  end while;
  s := F(s,myt, shortTableOfSigns);  Append(~res, [s]);
  return res;
end function;


function AGMPrime(b,t, shortTableOfSigns)
  R := Finfty(b,t, shortTableOfSigns);
  return R[#R][1][1][1];
end function;

/*
prec := 400;
C := ComplexField(prec);
val := C!0.866;
tau := Matrix(3,3, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, C!0.05+C.1*val/10, C!0.15+C.1*val/9, -C!0.4+C.1*val, C!0.2+C.1*val/20, C!0.05+C.1*val/10, C!0.2+C.1*val/20, C!0.1+C.1*val*C!1.5]);
thconstants := Matrix(8, 1, [ ThetaGenus3(i,tau)^2 : i in [0..7]]);
Finfty(thconstants,tau, ComputeOneTable(tau));
time AGMPrime(thconstants,tau, ComputeOneTable(tau));
*/


function AllDuplication(a)
  // Given theta_{0,b}(0,t) compute theta_{a,b}(0,t)

  // we still need to find a more genus g way to do this (Makarov = sum ab, not sum (-1)ab ab)
  n := #a;
  b := a;

  ThetaProducts := Matrix(8,8,
	[a[1]*b[1],a[1]*b[2],a[1]*b[3],a[1]*b[4],a[1]*b[5],a[1]*b[6],a[1]*b[7],a[1]*b[8],
	 a[2]*b[2],a[2]*b[1],a[2]*b[4],a[2]*b[3],a[2]*b[6],a[2]*b[5],a[2]*b[8],a[2]*b[7],
	 a[3]*b[3],a[3]*b[4],a[3]*b[1],a[3]*b[2],a[3]*b[7],a[3]*b[8],a[3]*b[5],a[3]*b[6],
	 a[4]*b[4],a[4]*b[3],a[4]*b[2],a[4]*b[1],a[4]*b[8],a[4]*b[7],a[4]*b[6],a[4]*b[5],
	 a[5]*b[5],a[5]*b[6],a[5]*b[7],a[5]*b[8],a[5]*b[1],a[5]*b[2],a[5]*b[3],a[5]*b[4],
	 a[6]*b[6],a[6]*b[5],a[6]*b[8],a[6]*b[7],a[6]*b[2],a[6]*b[1],a[6]*b[4],a[6]*b[3],
	 a[7]*b[7],a[7]*b[8],a[7]*b[5],a[7]*b[6],a[7]*b[3],a[7]*b[4],a[7]*b[1],a[7]*b[2],
	 a[8]*b[8],a[8]*b[7],a[8]*b[6],a[8]*b[5],a[8]*b[4],a[8]*b[3],a[8]*b[2],a[8]*b[1]]
		);
  hadam := HadamardMatrix(Parent(a[1]), 3);

  ThetaProducts := hadam*ThetaProducts;
  ThetaProducts := ElementToSequence(ThetaProducts/8);

  return ThetaProducts;
end function;

/*
prec := 60;
C := ComplexField(prec);
val := C!0.866;
tau := Matrix(3,3, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, C!0.05+C.1*val/10, C!0.15+C.1*val/9, -C!0.4+C.1*val, C!0.2+C.1*val/20, C!0.05+C.1*val/10, C!0.2+C.1*val/20, C!0.1+C.1*val*C!1.5]);
thconstants := [ ThetaGenus3(i,tau) : i in [0..7]];
time allthetas2tau := AllDuplication(thconstants);
for i := 1 to 64 do
  ComplexField(10)!Abs(allthetas2tau[i]-ThetaGenus3(i-1,2*tau)^2);
end for;
*/


function toInverse(b,t, tableOfSigns)
  // Given theta_i/theta_0, compute tau, and the equation linking the theta-constants or an extra equation that's linearly independent (so that it's a 7->7)
  // ya ya, it's weird to have a function which goal is to compute tau but you give it tau (for the signs) in the args
  p := Precision(Parent(b[1][1]));
  n := #Rows(b);
  CC := ComplexField(p);

  // First step: compute 1/theta00(z)^2, 1/theta00(0)^2
  theta000 := AGMPrime(b,t, tableOfSigns[1]);
  theta000 := 1/theta000;
  // then compute the other ones
  theta0withAequals0 := [];
  for i:=1 to n do Append(~theta0withAequals0, b[i][1]*theta000); end for;

  //theta0withAequals0;

  // Then compute everything at 2tau (simpler conceptually and generalizable to genus g)
  rootB := [ Sqrt(theta0withAequals0[i])*tableOfSigns[1][i] : i in [1..n]];
  sixtyfourThetaConstants := AllDuplication(rootB);
  
  // then give it to borchardt
  tt := 2*t;

  i3 := IdentityMatrix(CC,3);
  delta11 := Matrix(3,3, [CC!1,CC!0,CC!0,CC!0,CC!0,CC!0,CC!0,CC!0,CC!0]);
  delta22 := Matrix(3,3, [CC!0,CC!0,CC!0,CC!0,CC!1,CC!0,CC!0,CC!0,CC!0]);
  delta33 := Matrix(3,3, [CC!0,CC!0,CC!0,CC!0,CC!0,CC!0,CC!0,CC!0,CC!1]);
  delta12 := Matrix(3,3, [CC!0,CC!1,CC!0,CC!1,CC!0,CC!0,CC!0,CC!0,CC!0]);
  delta13 := Matrix(3,3, [CC!0,CC!0,CC!1,CC!0,CC!0,CC!0,CC!1,CC!0,CC!0]);
  delta23 := Matrix(3,3, [CC!0,CC!0,CC!0,CC!0,CC!0,CC!1,CC!0,CC!1,CC!0]);
  jm1 := (-tt-delta11)*(delta11*tt+delta11-i3)^(-1);
  jm2 := (-tt-delta22)*(delta22*tt+delta22-i3)^(-1);
  jm3 := (-tt-delta33)*(delta33*tt+delta33-i3)^(-1);
  jm12 := (-tt-delta12)*(delta12*tt+delta11+delta22-i3)^(-1);
  jm13 := (-tt-delta13)*(delta13*tt+delta11+delta33-i3)^(-1);
  jm23 := (-tt-delta23)*(delta23*tt+delta33+delta22-i3)^(-1);

  /*
  ComplexField(10)!Abs(ThetaGenus3(0,jm1)^2+C.1*tt[1][1]*ThetaGenus3(32,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(4,jm1)^2+C.1*tt[1][1]*ThetaGenus3(0,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(0,jm2)^2+C.1*tt[2][2]*ThetaGenus3(16,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(4,jm2)^2+C.1*tt[2][2]*ThetaGenus3(20,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(0,jm3)^2+C.1*tt[3][3]*ThetaGenus3(8,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(2,jm3)^2+C.1*tt[3][3]*ThetaGenus3(10,tt)^2);
  
  ComplexField(10)!Abs(ThetaGenus3(2,jm12)^2-(tt[1][2]^2-tt[1][1]*tt[2][2])*ThetaGenus3(32,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(6,jm12)^2-(tt[1][2]^2-tt[1][1]*tt[2][2])*ThetaGenus3(48,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(3,jm13)^2-(tt[1][3]^2-tt[1][1]*tt[3][3])*ThetaGenus3(34,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(5,jm13)^2-(tt[1][3]^2-tt[1][1]*tt[3][3])*ThetaGenus3(40,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(0,jm23)^2-(tt[2][3]^2-tt[2][2]*tt[3][3])*ThetaGenus3(0,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(3,jm23)^2-(tt[2][3]^2-tt[2][2]*tt[3][3])*ThetaGenus3(24,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(7,jm23)^2-(tt[2][3]^2-tt[2][2]*tt[3][3])*ThetaGenus3(28,tt)^2);
  */


  quotient8sur0 := sixtyfourThetaConstants[8+1]/sixtyfourThetaConstants[0+1];
  quotient16sur0 := sixtyfourThetaConstants[16+1]/sixtyfourThetaConstants[0+1];
  quotient24sur0 := sixtyfourThetaConstants[24+1]/sixtyfourThetaConstants[0+1];
  quotient32sur0 := sixtyfourThetaConstants[32+1]/sixtyfourThetaConstants[0+1];
  quotient40sur0 := sixtyfourThetaConstants[40+1]/sixtyfourThetaConstants[0+1];
  quotient48sur0 := sixtyfourThetaConstants[48+1]/sixtyfourThetaConstants[0+1];

  mu1 := AGMPrime(Matrix(8,1, [CC!1, sixtyfourThetaConstants[33+1]/sixtyfourThetaConstants[32+1], sixtyfourThetaConstants[34+1]/sixtyfourThetaConstants[32+1], sixtyfourThetaConstants[35+1]/sixtyfourThetaConstants[32+1], 1/quotient32sur0, sixtyfourThetaConstants[1+1]/sixtyfourThetaConstants[32+1], sixtyfourThetaConstants[2+1]/sixtyfourThetaConstants[32+1], sixtyfourThetaConstants[3+1]/sixtyfourThetaConstants[32+1]]), jm1, tableOfSigns[2]);
  mu2 := AGMPrime(Matrix(8,1, [CC!1, sixtyfourThetaConstants[17+1]/sixtyfourThetaConstants[16+1], 1/quotient16sur0, sixtyfourThetaConstants[1+1]/sixtyfourThetaConstants[16+1], sixtyfourThetaConstants[20+1]/sixtyfourThetaConstants[16+1], sixtyfourThetaConstants[21+1]/sixtyfourThetaConstants[16+1], sixtyfourThetaConstants[4+1]/sixtyfourThetaConstants[16+1], sixtyfourThetaConstants[5+1]/sixtyfourThetaConstants[16+1]]), jm2, tableOfSigns[3]);
  mu3 := AGMPrime(Matrix(8,1, [CC!1, 1/quotient8sur0, sixtyfourThetaConstants[10+1]/sixtyfourThetaConstants[8+1], sixtyfourThetaConstants[2+1]/sixtyfourThetaConstants[8+1], sixtyfourThetaConstants[12+1]/sixtyfourThetaConstants[8+1], sixtyfourThetaConstants[4+1]/sixtyfourThetaConstants[8+1], sixtyfourThetaConstants[14+1]/sixtyfourThetaConstants[8+1], sixtyfourThetaConstants[6+1]/sixtyfourThetaConstants[8+1]]), jm3, tableOfSigns[4]);

  mu12 := AGMPrime(Matrix(8,1, [CC!1, sixtyfourThetaConstants[1+1]/sixtyfourThetaConstants[0+1], quotient32sur0, sixtyfourThetaConstants[33+1]/sixtyfourThetaConstants[0+1], quotient16sur0, sixtyfourThetaConstants[17+1]/sixtyfourThetaConstants[0+1], quotient48sur0, sixtyfourThetaConstants[49+1]/sixtyfourThetaConstants[0+1]]), jm12, tableOfSigns[5]);
  mu13 := AGMPrime(Matrix(8,1, [CC!1, quotient32sur0, sixtyfourThetaConstants[2+1]/sixtyfourThetaConstants[0+1], sixtyfourThetaConstants[34+1]/sixtyfourThetaConstants[0+1], quotient8sur0, quotient40sur0, sixtyfourThetaConstants[10+1]/sixtyfourThetaConstants[0+1], sixtyfourThetaConstants[42+1]/sixtyfourThetaConstants[0+1]]), jm13, tableOfSigns[6]);
  mu23 := AGMPrime(Matrix(8,1, [CC!1, quotient16sur0, quotient8sur0, quotient24sur0, sixtyfourThetaConstants[4+1]/sixtyfourThetaConstants[0+1], sixtyfourThetaConstants[20+1]/sixtyfourThetaConstants[0+1], sixtyfourThetaConstants[12+1]/sixtyfourThetaConstants[0+1], sixtyfourThetaConstants[28+1]/sixtyfourThetaConstants[0+1]]), jm23, tableOfSigns[7]);

  // then extract the stuff
  computedtau11 := CC.1/(mu1*sixtyfourThetaConstants[32+1]);
  computedtau22 := CC.1/(mu2*sixtyfourThetaConstants[16+1]);
  computedtau33 := CC.1/(mu3*sixtyfourThetaConstants[8+1]);
  
  computedtau12 := 1/(mu12*sixtyfourThetaConstants[0+1]);
  computedtau13 := 1/(mu13*sixtyfourThetaConstants[0+1]);
  computedtau23 := 1/(mu23*sixtyfourThetaConstants[0+1]);
  
  //computedtau3 := Sqrt(computedtau3 + computedtau1*computedtau2);	// imaginary part positive (cause reduction)

  // because we used duplication formulas 
  computedtau11 := computedtau11/2; computedtau22 := computedtau22/2; computedtau33 := computedtau33/2;
  computedtau12 := computedtau12/4; computedtau13 := computedtau13/4; computedtau23 := computedtau23/4;

  // (We don't do this but maybe we should) compute the extra equation
  //eqBetweenGenus3ThetaConstants := CC!0;

  // compute a seventh quantity
  mydet := AGMPrime(Matrix(8,1, [CC!1, quotient8sur0, quotient16sur0, sixtyfourThetaConstants[24+1]/sixtyfourThetaConstants[0+1], quotient32sur0, quotient40sur0, quotient48sur0, sixtyfourThetaConstants[56+1]/sixtyfourThetaConstants[0+1]]), -tt^(-1), tableOfSigns[8]);
  mydet := -CC.1/(mydet*sixtyfourThetaConstants[0+1]); mydet := mydet/8;

  return [computedtau11, computedtau22, computedtau33, computedtau12, computedtau13, computedtau23, mydet];
end function;

/*
// check if we have the right action
prec := 150;
C := ComplexField(prec);
val := C!0.866;
tt := Matrix(3,3, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, C!0.05+C.1*val/10, C!0.15+C.1*val/9, -C!0.4+C.1*val, C!0.2+C.1*val/20, C!0.05+C.1*val/10, C!0.2+C.1*val/20, C!0.1+C.1*val*C!1.5]);

  i3 := IdentityMatrix(C,3);
  delta11 := Matrix(3,3, [C!1,C!0,C!0,C!0,C!0,C!0,C!0,C!0,C!0]);
  delta22 := Matrix(3,3, [C!0,C!0,C!0,C!0,C!1,C!0,C!0,C!0,C!0]);
  delta33 := Matrix(3,3, [C!0,C!0,C!0,C!0,C!0,C!0,C!0,C!0,C!1]);
  delta12 := Matrix(3,3, [C!0,C!1,C!0,C!1,C!0,C!0,C!0,C!0,C!0]);
  delta13 := Matrix(3,3, [C!0,C!0,C!1,C!0,C!0,C!0,C!1,C!0,C!0]);
  delta23 := Matrix(3,3, [C!0,C!0,C!0,C!0,C!0,C!1,C!0,C!1,C!0]);
  jm1 := (-tt-delta11)*(delta11*tt+delta11-i3)^(-1);
  jm2 := (-tt-delta22)*(delta22*tt+delta22-i3)^(-1);
  jm3 := (-tt-delta33)*(delta33*tt+delta33-i3)^(-1);
  jm12 := (-tt-delta12)*(delta12*tt+delta11+delta22-i3)^(-1);
  jm13 := (-tt-delta13)*(delta13*tt+delta11+delta33-i3)^(-1);
  jm23 := (-tt-delta23)*(delta23*tt+delta33+delta22-i3)^(-1);


  ComplexField(10)!Abs(ThetaGenus3(0,jm1)^2+C.1*tt[1][1]*ThetaGenus3(32,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(4,jm1)^2+C.1*tt[1][1]*ThetaGenus3(0,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(0,jm2)^2+C.1*tt[2][2]*ThetaGenus3(16,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(4,jm2)^2+C.1*tt[2][2]*ThetaGenus3(20,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(0,jm3)^2+C.1*tt[3][3]*ThetaGenus3(8,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(2,jm3)^2+C.1*tt[3][3]*ThetaGenus3(10,tt)^2);
  
  ComplexField(10)!Abs(ThetaGenus3(2,jm12)^2-(tt[1][2]^2-tt[1][1]*tt[2][2])*ThetaGenus3(32,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(6,jm12)^2-(tt[1][2]^2-tt[1][1]*tt[2][2])*ThetaGenus3(48,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(3,jm13)^2-(tt[1][3]^2-tt[1][1]*tt[3][3])*ThetaGenus3(34,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(5,jm13)^2-(tt[1][3]^2-tt[1][1]*tt[3][3])*ThetaGenus3(40,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(0,jm23)^2-(tt[2][3]^2-tt[2][2]*tt[3][3])*ThetaGenus3(0,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(3,jm23)^2-(tt[2][3]^2-tt[2][2]*tt[3][3])*ThetaGenus3(24,tt)^2);
  ComplexField(10)!Abs(ThetaGenus3(7,jm23)^2-(tt[2][3]^2-tt[2][2]*tt[3][3])*ThetaGenus3(28,tt)^2);


  ComplexField(10)!Abs(ThetaGenus3(2,-tt^(-1))^2-C.1*Determinant(tt)*ThetaGenus3(16,tt)^2);

*/


/*
prec := 400;
C := ComplexField(prec);
val := C!0.866;
tau := Matrix(3,3, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, C!0.05+C.1*val/10, C!0.15+C.1*val/9, -C!0.4+C.1*val, C!0.2+C.1*val/20, C!0.05+C.1*val/10, C!0.2+C.1*val/20, C!0.1+C.1*val*C!1.5]);
b := Matrix(8,1, [ThetaGenus3(0,tau)^2, ThetaGenus3(1,tau)^2, ThetaGenus3(2,tau)^2, ThetaGenus3(3,tau)^2, ThetaGenus3(4,tau)^2, ThetaGenus3(5,tau)^2, ThetaGenus3(6,tau)^2, ThetaGenus3(7,tau)^2]);
quotients0 := Matrix(8,1, [ b[i][1]/b[1][1] : i in [1..8]] );

time myres := toInverse(quotients0, tau, ComputeTableOfSigns(tau));

resultIwant := [ tau[1][1], tau[2][2], tau[3][3], tau[1][2]^2-tau[1][1]*tau[2][2], tau[1][3]^2-tau[1][1]*tau[3][3], tau[2][3]^2-tau[2][2]*tau[3][3], Determinant(tau)];

[ RealField(10)!Abs(resultIwant[i]-myres[i]) : i in [1..7]];
*/




function DiffFiniesOneStep(b, t,  tableOfSigns)
  // a is, as always, of size 4, but the first element is 1
  // Give approximations of the quotients (6 of them + 2 ones) with precision p and get an approximation with precision 2p
  prec := Precision(Parent(b[1][1]));
  CC := ComplexField(2*prec);

  newb := [ CC!b[i][1] : i in [1..#Rows(b)]];
  myBofSize7 := [ newb[i]/newb[1] : i in [2..#newb]];

  epsilon := CC!(10^(-prec));
  // we're going to assume that z and t are still good enough to determine the correct sign
  toInverseDeBase := toInverse( Matrix(8,1,[CC!1,myBofSize7[1],myBofSize7[2],myBofSize7[3],myBofSize7[4],myBofSize7[5],myBofSize7[6],myBofSize7[7]]) ,t, tableOfSigns);

  perturb := [];
  Append(~perturb, toInverse( Matrix(8,1,[CC!1,myBofSize7[1]+epsilon,myBofSize7[2],myBofSize7[3],myBofSize7[4],myBofSize7[5],myBofSize7[6],myBofSize7[7]]) ,t, tableOfSigns) );
  Append(~perturb, toInverse( Matrix(8,1,[CC!1,myBofSize7[1],myBofSize7[2]+epsilon,myBofSize7[3],myBofSize7[4],myBofSize7[5],myBofSize7[6],myBofSize7[7]]) ,t, tableOfSigns) );
  Append(~perturb, toInverse( Matrix(8,1,[CC!1,myBofSize7[1],myBofSize7[2],myBofSize7[3]+epsilon,myBofSize7[4],myBofSize7[5],myBofSize7[6],myBofSize7[7]]) ,t, tableOfSigns) );
  Append(~perturb, toInverse( Matrix(8,1,[CC!1,myBofSize7[1],myBofSize7[2],myBofSize7[3],myBofSize7[4]+epsilon,myBofSize7[5],myBofSize7[6],myBofSize7[7]]) ,t, tableOfSigns) );
  Append(~perturb, toInverse( Matrix(8,1,[CC!1,myBofSize7[1],myBofSize7[2],myBofSize7[3],myBofSize7[4],myBofSize7[5]+epsilon,myBofSize7[6],myBofSize7[7]]) ,t, tableOfSigns) );
  Append(~perturb, toInverse( Matrix(8,1,[CC!1,myBofSize7[1],myBofSize7[2],myBofSize7[3],myBofSize7[4],myBofSize7[5],myBofSize7[6]+epsilon,myBofSize7[7]]) ,t, tableOfSigns) );
  Append(~perturb, toInverse( Matrix(8,1,[CC!1,myBofSize7[1],myBofSize7[2],myBofSize7[3],myBofSize7[4],myBofSize7[5],myBofSize7[6],myBofSize7[7]+epsilon]) ,t, tableOfSigns) );

  // need to get an explicit form for the jacobian? bah, inverse(jacobian) is fine)
  jacobienne := ZeroMatrix(CC,7,7);
  for i:=1 to 7 do
    for j := 1 to 7 do
	jacobienne[i][j] := (perturb[j][i]-toInverseDeBase[i])/epsilon;
    end for;
  end for;
  //dispjac := jacobienne; dispjac := Matrix(6,6, [ ComplexField(3)!Abs(dispjac[i][j]) : i,j in [1..6]]); dispjac;

  changement := Matrix([ [ toInverseDeBase[1]-t[1][1], toInverseDeBase[2]-t[2][2], toInverseDeBase[3]-t[3][3], toInverseDeBase[4]-(t[1][2]^2-t[1][1]*t[2][2]), toInverseDeBase[5]-(t[1][3]^2-t[1][1]*t[3][3]), toInverseDeBase[6]-(t[2][3]^2-t[2][2]*t[3][3]), toInverseDeBase[7]-Determinant(t)] ]);

//  changement := Matrix([[ toInverseDeBase[1]-lambdaIwant[1], toInverseDeBase[2]-lambdaIwant[2], toInverseDeBase[3]-lambdaIwant[3],
//			  toInverseDeBase[4]-t[1][1], toInverseDeBase[5]-t[2][2], toInverseDeBase[6]-(t[1][2]^2-t[1][1]*t[2][2])]]);
  //RealField(10)!( Abs(changement[1][1]) + Abs(changement[1][2]) + Abs(changement[1][3]) + Abs(changement[1][4]) + Abs(changement[1][5]) + Abs(changement[1][6]) );
//  if ( RealField(10)!( Abs(changement[1][1]) + Abs(changement[1][2]) + Abs(changement[1][3]) + Abs(changement[1][4]) + Abs(changement[1][5]) + Abs(changement[1][6]) ) gt 10*Real(epsilon) ) then print "Something went wrong in the Newton"; end if;
  changement := changement*(Transpose(jacobienne)^(-1));
  //ComplexField(10)!( Abs(changement[1][1]) + Abs(changement[1][2]) + Abs(changement[1][3]) + Abs(changement[1][4]) + Abs(changement[1][5]) + Abs(changement[1][6]) );

  return Matrix(8,1,[1, myBofSize7[1]-changement[1][1], myBofSize7[2]-changement[1][2], myBofSize7[3]-changement[1][3], myBofSize7[4]-changement[1][4], myBofSize7[5]-changement[1][5], myBofSize7[6]-changement[1][6], myBofSize7[7]-changement[1][7]]);
end function;



/*
prec := 100;
C := ComplexField(2*prec);
CC:= ComplexField(prec);
val := C!0.866;
tau := Matrix(3,3, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, C!0.05+C.1*val/10, C!0.15+C.1*val/9, -C!0.4+C.1*val, C!0.2+C.1*val/20, C!0.05+C.1*val/10, C!0.2+C.1*val/20, C!0.1+C.1*val*C!1.5]);
val := CC!val;
tautau := Matrix(3,3, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, C!0.05+C.1*val/10, C!0.15+C.1*val/9, -C!0.4+C.1*val, C!0.2+C.1*val/20, C!0.05+C.1*val/10, C!0.2+C.1*val/20, C!0.1+C.1*val*C!1.5]);

//b := Matrix(8,1, [ThetaGenus3(0,tau)^2, ThetaGenus3(1,tau)^2, ThetaGenus3(2,tau)^2, ThetaGenus3(3,tau)^2, ThetaGenus3(4,tau)^2, ThetaGenus3(5,tau)^2, ThetaGenus3(6,tau)^2, ThetaGenus3(7,tau)^2]);
b := Matrix(8,1, [ThetaGenus3(0,tautau)^2, ThetaGenus3(1,tautau)^2, ThetaGenus3(2,tautau)^2, ThetaGenus3(3,tautau)^2, ThetaGenus3(4,tautau)^2, ThetaGenus3(5,tautau)^2, ThetaGenus3(6,tautau)^2, ThetaGenus3(7,tautau)^2]);


quotients0 := Matrix(8,1, [ b[i][1]/b[1][1] : i in [1..8]] );
//rrrrr := toInverse(quotients0, tau);
// [ rrrrr[1]-tau[1][1], rrrrr[2]-tau[2][2], rrrrr[3]-tau[3][3], rrrrr[4]-tau[1][2], rrrrr[5]-tau[1][3], rrrrr[6]-tau[2][3], rrrrr[7]];

apresdiffB := DiffFiniesOneStep(quotients0, tau, ComputeTableOfSigns(tau));



ComplexField(10)!Abs(apresdiffB[2][1]-ThetaGenus3(1,tau)^2/ThetaGenus3(0,tau)^2);
//C!quotients0[2][1]-ThetaGenus3(1,tau)^2/ThetaGenus3(0,tau)^2;

*/


LOW_PRECISION := 450;

function CalculThetas(tau)

  // Returns [theta_{0,b}(0,tau)]_{b=0..7}

  CC := Parent(tau[1][1]);

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
  taulow := Matrix(3,3, [Clow!tau[1][1], Clow!tau[1][2], Clow!tau[1][3], Clow!tau[2][1], Clow!tau[2][2], Clow!tau[2][3], Clow!tau[3][1],Clow!tau[3][2],Clow!tau[3][3]]);
  initB := NaiveThetaConstantsGenus3(taulow,1);
  if (flag eq 1) then
	return initB;
  else
	// we keep going
	initB := [ initB[i]^2 : i in [1..#initB]];
  end if;


  // we're going to need the signs of the quotients ; let's avoid recomputing them
  tableOfSigns := ComputeTableOfSigns(tau);




  b := Matrix(8,1, [ initB[i]/initB[1] : i in [1..8]]);
  // Finite differences to refine the result
  p := lowprec;
  while (p lt Precision(tau[1][1])) do
    p := 2*p;
    Cp := ComplexField(p);
    b := DiffFiniesOneStep(b, Matrix(3,3, [Cp!tau[1][1], Cp!tau[1][2], Cp!tau[1][3], Cp!tau[2][1], Cp!tau[2][2], Cp!tau[2][3], Cp!tau[3][1],Cp!tau[3][2],Cp!tau[3][3]]), tableOfSigns);
  end while;

  // We have the quotient of thetas, apply F to unstick them
  theta000 := AGMPrime(b,tau, tableOfSigns[1]);
  //theta00z, theta000;
  b := [ Sqrt(b[i][1]/theta000)*tableOfSigns[1][i] : i in [1..8]];

  b := [CC!b[1], CC!b[2], CC!b[3], CC!b[4], CC!b[5], CC!b[6], CC!b[7], CC!b[8]];
  return b;
end function;


function Magma8ThetaConstants(tau)
  theC := ComplexField(Precision(tau[1][1]));
  resB := [ThetaGenus3(0,tau), ThetaGenus3(1,tau), ThetaGenus3(2,tau), ThetaGenus3(3,tau), ThetaGenus3(4,tau), ThetaGenus3(5,tau), ThetaGenus3(6,tau), ThetaGenus3(7,tau)];
  return resB;
end function;


/*prec := 940;
C := ComplexField(prec);
val := C!0.866;
tau := Matrix(3,3, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, C!0.05+C.1*val/10, C!0.15+C.1*val/9, -C!0.4+C.1*val, C!0.2+C.1*val/20, C!0.05+C.1*val/10, C!0.2+C.1*val/20, C!0.1+C.1*val*C!1.5]);

time b := CalculThetas(tau);
time naiveB := NaiveThetaConstantsGenus3(tau,C!1);
//time expectedB := Magma8ThetaConstants(tau);
time notNeeded := ThetaGenus3(0,tau);

[ RealField(10)!Abs(b[i]-naiveB[i]) : i in [1..8]];
*/

/*
for i:= 0 to 10 do
	prec := 1000*2^i;
        printf "\nPrec : %o\n",prec;
	C := ComplexField(prec);
	val := C!0.866;
	tau := Matrix(3,3, [C!0.2+val*C.1/2, C!0.15+C.1*val/9, C!0.05+C.1*val/10, C!0.15+C.1*val/9, -C!0.4+C.1*val, C!0.2+C.1*val/20, C!0.05+C.1*val/10, C!0.2+C.1*val/20, C!0.1+C.1*val*C!1.5]);

	time b := CalculThetas(tau);
	time naiveB := NaiveThetaConstantsGenus3(tau,C!1);
	//time expectedA, expectedB := Magma8Thetas(tau);
	time notNeeded := ThetaGenus3(0,tau);
end for;

*/










