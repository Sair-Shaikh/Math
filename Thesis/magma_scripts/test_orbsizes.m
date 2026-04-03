import "magma_scripts/Utils.m" : ConvertToPolys;




F := Open("Dataset/orbits_point_counts.m", "r");
orbits := ReadObject(F);
delete F;

// Setup
F2 := GF(2);
G := GL(5, F2);
R<[x]>   := PolynomialRing(F2, 5);
V2, Bit2 := GModule(G, R, 2);
V3, Bit3 := GModule(G, R, 3);


for k -> orbit in orbits do 
    osize := orbit["size"];

    // Convert from int to polynomials in R<x>
    f2_int := orbit["f2int"];
    f3_int := orbit["f3int"];
    f2, f3 := ConvertToPolys(f2_int, f3_int, R, G, Bit2, Bit3);  // Now f2, f3 are in R
    f2, f3, osize;


    break;
end for;