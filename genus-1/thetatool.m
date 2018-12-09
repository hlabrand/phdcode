
/* No, really, I shouldn't have to write this. */
function parserealnumber(s, CC)
    foo:=s;
    k:=1;
    while k le #foo and foo[k] eq " " do k+:=1; end while;
    if k le #foo then
        foo:=Substring(foo, k, #foo+1-k); k:=1;
    else
        foo:="";
    end if;
    s:=1;
    if Regexp("^@NaN@", foo) then 
        error "output string contains Nan, unsupported: " cat foo;
    end if;
    t, seg, parts:=Regexp("^(-?[0-9]+)\\.([0-9]*)", foo);
    if t then
        a,b:=Explode(parts);
        m:=StringToInteger(a);
        s:=a[1] eq "-" select -1 else 1;        // avoid 0
        m*:=s;
        m+:=StringToInteger(b)/10^#b;
    else
        t, seg, parts:=Regexp("^(-?[0-9]+)", foo);
        if t then
            a:=Explode(parts);
            m:=StringToInteger(a);
            s:=a[1] eq "-" select -1 else 1;        // avoid 0
            m*:=s;
        else
            t, seg, parts:=Regexp("(-?)\\.([0-9]*)", foo);
            assert t; /* we're in trouble here */
            a,b:=Explode(parts);
            m:=StringToInteger(b)/10^#b;
            s:=StringToInteger(a cat "1");
        end if;
    end if;
    m*:=s;
    k+:=#seg;
    if k le #foo then
        foo:=Substring(foo, k, #foo+1-k); k:=1;
    else
        foo:="";
    end if;
    t, seg, parts := Regexp("^[eE](-?[0-9]+)", foo);
    if t then
        e:=StringToInteger(parts[1]);
        m*:=10^e;
        k+:=#seg;
        if k le #foo then
            foo:=Substring(foo, k, #foo+1-k); k:=1;
        else
            foo:="";
        end if;
    end if;
    return CC!m, foo;
end function;

function output_to_complex_vector(output, CC)
    k:=1;
    L:=[];
    foo:=output;
    tm:=0;
    while k le #foo do
        if foo[k] ne "(" then
            /* timing info */
            foo:=Substring(foo, k, #foo+1-k); k:=1;
            tm, foo:=parserealnumber(foo, CC);
            assert foo[k] eq "\n"; k+:=1;
            assert k gt #foo;
            break;
        end if;
        assert foo[k] eq "(";
        k+:=1;
        foo:=Substring(foo, k, #foo+1-k); k:=1;
        rr, foo:=parserealnumber(foo, CC);
        ii, foo:=parserealnumber(foo, CC);
        assert foo[k] eq ")"; k+:=1;
        assert foo[k] eq "\n"; k+:=1;
        Append(~L, CC![rr,ii]);
    end while;
    return Vector(L), tm;
end function;

function input_from_complex_vector(vec)
    return &cat [Sprintf("(%o %o)\n", Real(t), Imaginary(t)) : t in Eltseq(vec)];
end function;





// These are the functions you're looking for.

intrinsic c_FastThetas(z :: FldComElt, tau :: FldComElt) -> SeqEnum
  { Fast computation of all 6 thetas. }
  CC := Parent(z);
  input := input_from_complex_vector([z, tau]);
  output := Pipe( "fastthetas " cat "time FastThetas" cat " " cat IntegerToString(BitPrecision(CC)), input);
  //output;
  res := output_to_complex_vector(output, CC);
  return res; 
end intrinsic;

intrinsic c_DupontThetaConstants(tau :: FldComElt) -> SeqEnum
  { Fast computation of all 3 theta-constants. }
  CC := Parent(tau);
  input := input_from_complex_vector([tau]);
  output := Pipe( "fastthetas " cat "time DupontThetaConstants" cat " " cat IntegerToString(BitPrecision(CC)), input);
  //output;
  res := output_to_complex_vector(output, CC);
  return res;
end intrinsic;

intrinsic c_NaiveTheta(z :: FldComElt, tau :: FldComElt) -> SeqEnum
  { Compute theta(z,tau). }
  CC := Parent(z);
  input := input_from_complex_vector([z, tau]);
  output := Pipe( "fastthetas " cat "time NaiveTheta" cat " " cat IntegerToString(BitPrecision(CC)), input);
  res := output_to_complex_vector(output, CC);
  return res;
end intrinsic;

intrinsic c_NaiveThetaFunctionAndThetaConstants_00_01(z :: FldComElt, tau :: FldComElt) -> SeqEnum
  { Compute theta[00,01]([z,0],tau). }
  CC := Parent(z);
  input := input_from_complex_vector([z, tau]);
  output := Pipe( "fastthetas " cat "time NaiveThetaFunctionAndThetaConstants_00_01" cat " " cat IntegerToString(BitPrecision(CC)), input);
  res := output_to_complex_vector(output, CC);
  return res;
end intrinsic;

intrinsic c_NaiveThetaFunctionAndThetaConstants_00_01_10(z :: FldComElt, tau :: FldComElt) -> SeqEnum
  { Compute theta[00,01,10]([z,0],tau). }
  CC := Parent(z);
  input := input_from_complex_vector([z, tau]);
  output := Pipe( "fastthetas " cat "time NaiveThetaFunctionAndThetaConstants_00_01_10" cat " " cat IntegerToString(BitPrecision(CC)), input);
  res := output_to_complex_vector(output, CC);
  return res;
end intrinsic;













// Other functions



intrinsic c_F(L :: SeqEnum) -> SeqEnum
  { Applies the function F to 4 numbers. }
  CC:=Parent(L[1]);
  input := input_from_complex_vector(L);
  output := Pipe( "fastthetas " cat "time AGMPrime" cat " " cat IntegerToString(BitPrecision(CC)), input);
  res := output_to_complex_vector(output, CC);
  return res;
end intrinsic;

intrinsic c_Pow2PowN(z :: FldComElt, n :: RngIntElt) -> FldComElt
  { Computes x^2^n.}
  CC:=Parent(z);
  input := input_from_complex_vector([z,CC!n]);
  output := Pipe( "fastthetas " cat "time Pow2PowN" cat " " cat IntegerToString(BitPrecision(CC)), input);
  res := output_to_complex_vector(output, CC);
  return res;
end intrinsic;

intrinsic c_Finfty4args(a :: SeqEnum) -> SeqEnum
  { Computes Finfty(a[1],a[2],a[3],a[4]) ; returns lambda and mu. }
  CC := Parent(a[1]);
  input := input_from_complex_vector(a);
  output := Pipe( "fastthetas " cat "time Finfty4args" cat " " cat IntegerToString(BitPrecision(CC)), input);
  res := output_to_complex_vector(output, CC);
  return res;
end intrinsic;

intrinsic c_Finfty2args(quo1 :: FldComElt, quo2 :: FldComElt) -> SeqEnum
  { Computes Finfty(1, quo1, 1, quo2) ; returns lambda and mu. }
  CC := Parent(quo1);
  input := input_from_complex_vector([quo1,quo2]);
  output := Pipe( "fastthetas " cat "time Finfty2args" cat " " cat IntegerToString(BitPrecision(CC)), input);
  //output;
  res := output_to_complex_vector(output, CC);
  return res;
end intrinsic;

intrinsic c_WholeFunction(quo1 :: FldComElt, quo2 :: FldComElt) -> SeqEnum
  { Computes the lambda and mu corresponding to the input theta01/theta00([0,z],tau). }
  CC := Parent(quo1);
  input := input_from_complex_vector([quo1,quo2]);
  output := Pipe( "fastthetas " cat "time WholeFunction" cat " " cat IntegerToString(BitPrecision(CC)), input);
  res := output_to_complex_vector(output, CC);
  return res;
end intrinsic;

intrinsic c_NewtonStep(quo1 :: FldComElt, quo2 :: FldComElt, lambda :: FldComElt, mu :: FldComElt) -> SeqEnum
  { One step of finite differences : if you give it approx of the quotients at prec p + lambda and mu at prec 2p, it gives you an approx of the quotients at prec 2p. }
  CC := Parent(quo1);
  input := input_from_complex_vector([quo1,quo2, lambda, mu]);
  output := Pipe( "fastthetas " cat "time NewtonStep" cat " " cat IntegerToString(BitPrecision(CC)), input);
  res := output_to_complex_vector(output, CC);
  return res;
end intrinsic;





intrinsic c_InternalNaiveThetaFunctionAndThetaConstants_00_01(z :: FldComElt, tau :: FldComElt) -> SeqEnum
  { Compute theta[00,01]([z,0],tau). }
  CC := Parent(z);
  input := input_from_complex_vector([z, tau]);
  output := Pipe( "fastthetas " cat "time InternalNaiveThetaFunctionAndThetaConstants_00_01" cat " " cat IntegerToString(BitPrecision(CC)), input);
  res := output_to_complex_vector(output, CC);
  return res;
end intrinsic;

intrinsic c_NaiveThetaConstants_00_01(tau :: FldComElt) -> SeqEnum
  { Compute theta[00,01](0,tau) with small precision (40 bits). }
  CC := Parent(tau);
  input := input_from_complex_vector([tau]);
  output := Pipe( "fastthetas " cat "time NaiveThetaConstants_00_01" cat " " cat IntegerToString(BitPrecision(CC)), input);
  res := output_to_complex_vector(output, CC);
  return res;
end intrinsic;

intrinsic c_NaiveThetaFunctionAndThetaConstantsSmallPrec_00_01(z :: FldComElt, tau :: FldComElt) -> SeqEnum
  { Compute theta[00,01]([z,0],tau) with small precision (40 bits). }
  return c_NaiveThetaFunctionAndThetaConstants_00_01(z, ComplexField(14)!tau);
end intrinsic;

intrinsic c_NaiveThetaSmallPrec(z :: FldComElt, tau :: FldComElt) -> SeqEnum
  { Compute theta with small precision (40 bits). }
  return c_NaiveTheta(z, ComplexField(14)!tau);
end intrinsic;

intrinsic c_NaiveThetaFunctionAndThetaConstantsSmallPrec_00_01_10(z :: FldComElt, tau :: FldComElt) -> SeqEnum
  { Compute theta[00,01,10]([z,0],tau) with small precision (40 bits). }
  return c_NaiveThetaFunctionAndThetaConstants_00_01_10(z, ComplexField(14)!tau);
end intrinsic;

intrinsic c_NaiveThetaConstantsSmallPrec_00_01(tau :: FldComElt) -> SeqEnum
  { Compute theta[00,01](0,tau) with small precision (40 bits). }
  return c_NaiveThetaConstants_00_01(ComplexField(14)!tau);
end intrinsic;

