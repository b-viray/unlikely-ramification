/******************************************************************************/
//Determine the relations on the polynomial given by coefs for that polynomial to be a scalar times a perfect power.
//coefs starts with the constant term; m >= 2
function GetMthPowerRelns(m, coefs)
    assert #coefs eq m+1;
    assert m ge 2;

    RelnMatrix := [];
    for i := 0 to m-1 do
        Append(~RelnMatrix, [Integers()!((m - i)/GCD(m-i, i+1))*coefs[i+1], Integers()!((i+1)/GCD(m-i, i+1))*coefs[i+2]]);
    end for;

    RelnMatrix := Matrix(RelnMatrix);//Reln Matrix is rank 1 and every row is a scalar multiple of [-alpha, 1], where alpha is the only root
    //so we may find the linear form unless every row is zero. For a nonzero coefs vector, if every row of RelnMatrix is zero then the base ring has characteristic dividing m.

    Relns := ideal<Parent(coefs[1]) | Minors(RelnMatrix, 2)>;

    return Relns;

end function;
/******************************************************************************/

//Takes as input: A genus 2 curve C, with an effective integral degree 3 divisor D3, such that 
// 1) 2*D3 = 3*K_C,  
// 2) D3 - K_C is ineffective,
// 3) The torsion packet of a Weierstrass point of C consists only of the Weierstrass points

//The linear systems |K_C|, |D_3| embed C as a (2,3) curve in P^1\times P^1. 

//The function computes the odd primes such that there exist hyperplanes in |nK_C + D_3| defined over F_p that are totally tangent to C at an F_p-point. For n>0, this can be viewed as intersecting with an irreducible (1,n) curve

function FpPointInYp(coefs_tuple, n, g, p)
    assert #coefs_tuple eq 3;
    assert [#c : c in coefs_tuple] eq [g+2, g+2, g+2];
    coefs0, coefs1, coefs2 := Explode(coefs_tuple);

    if n eq 0 then
        A<T> := PolynomialRing(GF(p),1);
        new_poly := [coefs0[i]*T^2 + coefs1[i]*T + coefs2[i] : i in [1..g+2]];
        I := GetMthPowerRelns(g + 1, new_poly);
        Z := Spec(quo<A | I>);
        ratpts := RationalPoints(Z);
        print p, ratpts;
        return #ratpts gt 0; 
    else
        A := PolynomialRing(GF(p), 2*(n+1));
        AssignNames(~A,  ["a" cat IntegerToString(i) : i in [0..n]] cat ["b" cat IntegerToString(i) : i in [0..n]]);
        as := [A.i : i in [1..n+1]];
        bs := [A.i : i in [n+2..2*(n+1)]];

        B<x> := PolynomialRing(A);

        h0 := Polynomial([A!c : c in coefs0]);
        h1 := Polynomial([A!c : c in coefs1]);
        h2 := Polynomial([A!c : c in coefs2]);
        a := Polynomial([A!c : c in as]);
        b := Polynomial([A!c : c in bs]);
        //need to ensure that as and bs have no common factor
        Res := Resultant(a,b);

        I := GetMthPowerRelns(2*n + g + 1, Coefficients(h0*b^2 +h1*a*b + h2*a^2));
        I := Saturation(I, ideal<A | Res>);
        Z := Proj(quo<A | I>);

        ratpts := RationalPoints(Z);
        print p, ratpts;
        return #ratpts gt 0; 
    end if;
end function;


function ComputeTotallyTamelyRamifiedPrimes_Affine(coefs_tuple, n, g)
    assert #coefs_tuple eq 3;
    assert [#c : c in coefs_tuple] eq [g+2, g+2, g+2];
    coefs0, coefs1, coefs2 := Explode(coefs_tuple);

    if n eq 0 then
        A<T> := PolynomialRing(Integers(),1);
        new_poly := [coefs0[i]*T^2 + coefs1[i]*T + coefs2[i] : i in [1..g+2]];
        I := GetMthPowerRelns(g + 1, new_poly);
        J := EliminationIdeal(I, {});
               TRPrimes := PrimeDivisors(Integers()!Generators(J)[1]);
        RealizedTRPrimes := [];
        for p in TRPrimes do
            if p ne 2 and FpPointInYp(coefs_tuple, n, g, p) then
                Append(~RealizedTRPrimes, p);
            end if;
        end for;
        
        return Sort(SetToSequence(SequenceToSet(RealizedTRPrimes)));
    else
        A := PolynomialRing(Integers(), 2*(n+1));
        AssignNames(~A,  ["a" cat IntegerToString(i) : i in [0..n]] cat ["b" cat IntegerToString(i) : i in [0..n]]);
        as := [A.i : i in [1..n+1]];
        bs := [A.i : i in [n+2..2*(n+1)]];

        B<x> := PolynomialRing(A);

        h0 := Polynomial([A!c : c in coefs0]);
        h1 := Polynomial([A!c : c in coefs1]);
        h2 := Polynomial([A!c : c in coefs2]);
        a := Polynomial([A!c : c in as]);
        b := Polynomial([A!c : c in bs]);
        //need to ensure that as and bs have no common factor
        Res := Resultant(a,b);

        I := GetMthPowerRelns(2*n + g + 1, Coefficients(h0*b^2 +h1*a*b + h2*a^2));
        I := Saturation(I, ideal<A | Res>);

        J := EliminationIdeal(I, {});
        TRPrimes := PrimeDivisors(Integers()!Generators(J)[1]);
        RealizedTRPrimes := [];
        // print "computed TRPrimes";
        for p in TRPrimes do
        // print "testing ", p;
            if p ne 2 and FpPointInYp(coefs_tuple, n, g, p) then
                Append(~RealizedTRPrimes, p);
            end if;
        end for;

        return Sort(SetToSequence(SequenceToSet(RealizedTRPrimes)));
    end if;
end function;