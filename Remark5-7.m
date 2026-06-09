R<z0,z1,z2> := PolynomialRing(Rationals(), 3);
A<t> := PolynomialRing(Rationals());

f := (z0^2 + z1^2)^2 + z0*z2*(z1^2 + z2^2);

X := Proj(quo<R|f>);

assert IsNonsingular(X);

//intersection with ell = z0 - 50*z2

F := Evaluate(f, [50, t, 1]);

K := NumberField(F);
OK := MaximalOrder(K);
five := ideal<OK | 5>;
factors := Factorization(five);
assert #factors eq 1;
assert factors[1][2] eq 2;

