R<x> := PolynomialRing(Rationals()); 
C := HyperellipticCurve(R![-1, -4, -2, 2, -6, 1, -1], R![0, 1, 1]); //LMFDB curve 7676
cond := 7676;
X,pi:= SimplifiedModel(C);
WD := Divisor(X, Scheme(X, X.2));
HC := Divisor(X, Scheme(X, X.1));
D := Divisor(Support(WD)[1]) - HC;

p := 31;
assert cond mod p ne 0;
Xp := ChangeRing(X, GF(p));
Jp := Jacobian(Xp);
XpPoints := RationalPoints(Xp);
HCp := Divisor(Xp, Scheme(Xp, Xp.1));
WDp :=Divisor(Xp, Scheme(Xp, Xp.1^3 + 5*Xp.1*Xp.3^2 + 2*Xp.3^3));
Dp := -HCp;
for pp in Support(WDp) do
    Dp +:= Divisor(pp);
end for;

XpWPs := [P : P in XpPoints | P[2] eq 0];
XpNonWP := [P : P in XpPoints | P[2] ne 0]; //points that are not Weierstrass points
N := #Jp;
odd := Integers()!(N/2^Valuation(N,2));
HasTotallyRamifiedOddExtn := false;


for P in XpNonWP do //looping over non-Weierstrass points
    if IsPrincipal(odd*(Dp - Divisor(P))) then //checking if D- p has odd order.
        HasTotallyRamifiedOddExtn := true;
    end if;
    assert IsPrincipal(N*(Dp - Divisor(P)));
end for;

assert not HasTotallyRamifiedOddExtn; //we find that there are no totally ramified odd extensions at 31, which is consistent with the fact that 31 is not in Td for this curve.
