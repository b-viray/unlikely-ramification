load "Genus2OddDegreeTotallyRamifiedPrimes.m";

GoodCurves := [
    //curves from the LMFDB that have trivial MW group, an effective of divisor of degree 3 and that the torsion packet consists only of Weierstrass points (using Bjorn Poonen's pari code)
<"7676.a.7676.1", 7676, [[-1, -4, -2, 2, -6, 1, -1], [0, 1, 1]]>,
<"9977.a.9977.1", 9977, [[-1, -5, -9, -2, -5, 0, -1], [1, 0, 1]]>,
<"11259.a.11259.1", 11259, [[-1, 1, -2, 5, -6, 7, -8], [1, 1]]>,
<"13072.a.13072.1", 13072, [[-1, -2, 0, 2, -1, -2, -1], [1]]>,
<"15504.a.62016.1", 15504, [[-1, 3, -5, 1, 6, -2, -3], [1]]>,
<"19170.c.76680.1", 19170, [[-68, -82, -30, 8, 1, -3, -1], [0, 1, 0, 1]]>,
<"21989.a.21989.1", 21989, [[-16, 33, -17, -8, 8, 0, -1], [0, 0, 1]]>,
<"22032.a.22032.1", 22032, [[-1, 0, 1, 2, 0, -2, -1], [1]]>,
<"22185.a.22185.1", 22185, [[-1, 1, 0, 1, -3, -3, -1], [1, 1, 0, 1]]>,
<"23033.a.23033.1", 23033, [[-1, 0, 1, 0, -1, -1, -1], [1, 1, 0, 1]]>,
<"23735.a.23735.1", 23735, [[-1, 12, -58, 46, 0, -9, 2], [1, 1]]>,
<"25707.a.77121.1", 25707, [[-12, 7, -6, -5, 1, -1, -1], [0, 1, 0, 1]]>,
<"29636.a.474176.1", 29636, [[-1, 3, -5, 0, -1, 2, -1], [1, 1, 0, 1]]>,
<"31174.a.342914.1", 31174, [[-6, -20, -4, 17, -14, 4, -1], [0, 1, 1]]>,
<"31668.a.63336.1", 31668, [[-3, 12, -18, 5, 4, -3, -1], [0, 1, 1]]>,
<"36114.a.36114.1", 36114, [[-1, 2, 0, -5, 5, -15, 12], [1, 1]]>,
<"41151.a.41151.1", 41151, [[-4, -17, -20, 1, 9, -1, -1], [1, 1, 0, 1]]>,
<"43335.b.43335.1", 43335, [[-1, 1, -2, -4, 4, 0, -1], [1, 1, 0, 1]]>,
<"45175.a.45175.1", 45175, [[-1, 0, 12, -28, 18, -9, 2], [1, 1]]>,
<"47449.a.47449.1", 47449, [[-1, 0, 1, -3, 0, 1, -1], [1, 1, 0, 1]]>
];

print "===========================================================";
print "Computing totally ramified primes in D3 and KC + D3 for several genus 2 curves";
print "============================================================";
RamificationOutput := [<"test", 29, [[1,0,0,0,0,0],[1,0,0]], [3], [5]>];

for Cdesc in GoodCurves do
    polys := Cdesc[3];
    cond := Cdesc[2];
    h := polys[2];
    g := polys[1];
    factors := Factorization(Polynomial(h)^2 + 4*Polynomial(g));
    assert #factors eq 2;

    h0 := factors[1][1];
    h2 := Parent(h0)!(-(Polynomial(h)^2 + 4*Polynomial(g))/h0);
    assert h0*(-h2) eq Polynomial(h)^2 + 4*Polynomial(g);

    print Cdesc[1];
    // print "Primes of bad reduction (and 2):", Sort(PrimeDivisors(2*Cdesc[2]));


    coefs_tuple := [Coefficients(h0), [Parent(h[1])!0, Parent(h[1])!0, Parent(h[1])!0, Parent(h[1])!0], Coefficients(h2)];
    TRprimes0 := ComputeTotallyTamelyRamifiedPrimes_Affine(coefs_tuple, 0, 2);
    // print "Totally ramified primes for D3", TRprimes0;

    TRprimes1 := ComputeTotallyTamelyRamifiedPrimes_Affine(coefs_tuple, 1, 2);
    // print "Totally ramified primes for KC + D3", TRprimes1;


    // TRprimes2 := ComputeTotallyTamelyRamifiedPrimes_Affine(coefs_tuple, 2, 2);// n = 2 times out
    print "-----------------------------------";
    L := <
        Cdesc[1], Cdesc[2], Cdesc[3], 
        [p : p in TRprimes0 | 2*Cdesc[2] mod p ne 0], 
        [p : p in TRprimes1 | 2*Cdesc[2] mod p ne 0]
        >;
    print L;
    Append(~RamificationOutput, L); 
    // print RamificationOutput;
end for;


print "=========================================================";
print RamificationOutput;
exit;
