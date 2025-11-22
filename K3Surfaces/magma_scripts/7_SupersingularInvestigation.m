
import "magma_scripts/Utils.m" : ConvertToPolys;

// Select supersingular indices
N := 20;
F := Open("Dataset/zeta_functions/zeta_fails.m", "r");
selected_orbits := ReadObject(F);
delete F;
selected_orbits := selected_orbits[N];

// Load dataset
F := Open("Dataset/orbits_point_counts.m", "r");
orbits := ReadObject(F);
delete F;


GetPointCountAndKnownFactor := function(orbit) 
    R<T> := PolynomialRing(Rationals());
    assert IsDefined(orbit, "pt_count");
    assert IsDefined(orbit, "contains_line");

    pt_count := orbit["pt_count"];

    if orbit["contains_line"] eq 1 then
        known_factor :=  (T-2)^2;
    else 
        known_factor := (T-2)^1;
    end if;

    return pt_count, known_factor;

end function;


HalfWeil := function(point_counts, known_factor) 

    R<T> := PolynomialRing(Rationals());
    assert #point_counts ge 11;

    traces := [(point_counts[i] - 1 - 2^(2*i)) : i in [1..#point_counts]];

    wps := FrobeniusTracesToWeilPolynomials(traces, 2, 2, 22 : KnownFactor := known_factor);

    if #wps eq 1 then return wps; end if;

    wps := [wp : wp in wps | CheckWeilPolynomial(wp, 2, 1 : SurfDeg := 6)];

    return wps;
end function;




// Setup
F2 := GF(2);
G := GL(5, F2);
R<[x]>       := PolynomialRing(F2, 5);
V2, Bit2 := GModule(G, R, 2);
V3, Bit3 := GModule(G, R, 3);

for k -> orbit in orbits do 
    if not k in selected_orbits then continue; end if;

    // Convert from int to polynomials in R<x>
    f2_int := orbit["f2int"];
    f3_int := orbit["f3int"];
    f2, f3 := ConvertToPolys(f2_int, f3_int, R, G, Bit2, Bit3);  // Now f2, f3 are in R

    pt_count, known_factor := GetPointCountAndKnownFactor(orbit);

    char_polys := HalfWeil(pt_count, known_factor);

    [WeilPolynomialToRankBound(cand, 2) : cand in char_polys];


end for;
