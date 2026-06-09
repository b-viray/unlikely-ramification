load "Genus2OddDegreeTotallyRamifiedPrimes.m";


// https://www.lmfdb.org/Genus2Curve/Q/52969/b/52969/1
//Genus 2 curve 52969.b.52969.1
//has 3 torsion point in Weierstrass torsion packet
    h := [Integers()!1, 1, 0, 1];
    g := [Integers()!-1, 0, -1, 0, -1, -2, -1];
    factors := Factorization(Polynomial(h)^2 + 4*Polynomial(g));
    assert #factors eq 2;

    h0 := factors[1][1];
    h2 := Parent(h0)!(-(Polynomial(h)^2 + 4*Polynomial(g))/h0);
    assert h0*(-h2) eq Polynomial(h)^2 + 4*Polynomial(g);

    coefs_tuple := [Coefficients(h0), [Parent(h[1])!0, Parent(h[1])!0, Parent(h[1])!0, Parent(h[1])!0], Coefficients(h2)];

    coefs0, coefs1, coefs2 := Explode(coefs_tuple);

    A<T> := PolynomialRing(Integers(),1);
    new_poly := [coefs0[i]*T^2 + coefs1[i]*T + coefs2[i] : i in [1..4]];
    I := GetMthPowerRelns(3, new_poly);
    assert (T^2 + 7)^2 in I; // have totally ramified degree 3 points at all primes where T^2 + 7 has a root.
    //This is because -14*h0 + 2*h2 is a perfect cube
