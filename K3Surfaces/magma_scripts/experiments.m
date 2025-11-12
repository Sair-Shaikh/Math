import "magma_scripts/Utils.m" : ConvertToPolys;

F := Open("Dataset/orbits_lines_secrtpt.m", "r");
orbits := ReadObject(F);
delete F;
"Smooth Surfaces With A Secant Line: ", #orbits;

// Setup
F2 := GF(2);
G := GL(5, F2);

R<[x]>       := PolynomialRing(F2, 5);
Fiber<a,b,c> := PolynomialRing(F2, 3);
R3<u,v,w>    := PolynomialRing(Fiber, 3);
R2<s, t>     := PolynomialRing(Fiber, 2);
RP1<y,z>     := PolynomialRing(F2, 2);

V2, Bit2 := GModule(G, R, 2);
V3, Bit3 := GModule(G, R, 3);

new_orbits := AssociativeArray();
count := 0;
for key in Keys(orbits) do

    orbit := orbits[key];

    // Convert from int to polynomials in R<x>
    f2_int := orbit["f2int"];
    f3_int := orbit["f3int"];
    f2, f3 := ConvertToPolys(f2_int, f3_int, R, G, Bit2, Bit3);  // Now f2, f3 are in R


    if f2_int eq 1281 then 
        perm := [ x[2], x[1], x[3], x[4], x[5]];
        f2 := Evaluate(f2, perm);
        f3 := Evaluate(f3, perm);
    end if;


    orbit["f2int"] := Seqint([Integers() !Bit2(f2)[i] : i in [1..#Basis(V2)]], 2);
    orbit["f3int"] := Seqint([Integers() !Bit3(f3)[i] : i in [1..#Basis(V3)]], 2);
    new_orbits[key] := orbit;

    count +:= 1;
    if count mod 10000 eq 0 then count; end if;

end for;

F := Open("Dataset/orbits_lines_secrtpt2.m", "w");
WriteObject(F, new_orbits);
delete F;
