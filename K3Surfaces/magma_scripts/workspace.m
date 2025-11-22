import "magma_scripts/Utils.m" : PolyToBitVector, PolyToCpp;

R<H,a,b,c> := PolynomialRing(GF(2), 4, "lex");

str_reps := Reverse([PolyToCpp(Evaluate(mon, [1, a, b, c]), 2) : mon in MonomialsOfDegree(R, 6)]);

str := "#define MONS \\\n";
SetColumns(0);
for i in [1..#str_reps] do 
    rep := str_reps[85-i];
    str cat:= Sprintf("mons[%o] = ", i-1) cat rep;
    if not i eq #str_reps then 
        str cat:= ",\\\n";
    else
        str cat:= ";";
    end if;
end for;


arr := [5, 10, 15];
mons := Reverse([mon : mon in MonomialsOfDegree(R, 6)]);

R2<a,b,c> := PolynomialRing(GF(2), 3);
f := a^3 + a^2*c + b^3;
f;
PolyToBitVector(f);
PolyToBitVector(1);
PolyToBitVector(0);


// fname := "CppLib/monomials_bigq.h";
// F := Open(fname, "w");
// delete F; 
// PrintFile(fname, str);
