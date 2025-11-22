
import "magma_scripts/Utils.m" : ConvertToPolys;
/* 
    We will use the Weil conjectures to convert the point counts into the Weil Polynomials. 

    The Hodge diamond for K3 surfaces gives, alongside the Rationality and Betti numbers conjecture: 
    \zeta(X, s) = P_2(T)/(1-T)(1-q^2T), T = q^-s.

    Thus, we only need to determine P_2(T).
    We know that P_2(T) is the characteristic polynomial of the Frobenius acting on H^2_et(X_bar, \Q_l). 

    By the Leftschetz Trace Formula, 
    #X(F_q^m) = \sum_{i=0}^4 (-1)^i Tr(Frob^m | H^i) = 1 + Tr(Frob^m | H^2) + q^2m
    => Tr(From^m | H^2) = #X(F_q^m) - 1 - q^2m

    Then you use Newton's Identities to convert these: 
    Tr(From^m | H^2) = \sum \alpha^m

    Using the functional equation, we can determine #X(F_q^n) for n = 12, ..., 22 using #X(F_q^n) for n = 1, ..., 11.
*/

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

    // if #wps eq 1 then return wps; end if;
    // // Check if rk(Pic) = 2 and in that case enforce multiplicity of q > multiplicity of -q
    // wps := [ wp : wp in wps | (not WeilPolynomialToRankBound(wp, 2) eq 2) or (Valuation(wp, T-2) gt Valuation(wp, T+2)) ];

    return wps;
end function;

F := Open("Dataset/orbits_point_counts.m", "r");
orbits := ReadObject(F);
delete F;

// Setup
F2 := GF(2);
G := GL(5, F2);
R<[x]>   := PolynomialRing(F2, 5);
V2, Bit2 := GModule(G, R, 2);
V3, Bit3 := GModule(G, R, 3);

prog := 0;
failures := AssociativeArray();
for key in Keys(orbits) do 
    orbit := orbits[key];

    // Convert from int to polynomials in R<x>
    f2_int := orbit["f2int"];
    f3_int := orbit["f3int"];
    f2, f3 := ConvertToPolys(f2_int, f3_int, R, G, Bit2, Bit3);

    pt_count, known_factor := GetPointCountAndKnownFactor(orbit);

    char_polys := HalfWeil(pt_count, known_factor);

    assert #char_polys gt 0;

    if #char_polys eq 2 then 
        // Find lowest index above 11 for which coeff is 0
        coeffs := Coefficients(char_polys[1]);
        // #coeffs, coeffs;
        assert coeffs[12] eq 0;

        for i in [13..22] do 
            if not coeffs[i] eq 0 then 
                if not IsDefined(failures, i-1) then 
                    failures[i-1] := [];
                end if;
                Append(~failures[i-1], key);
                break;
            end if;
            if not IsDefined(failures, i-1) then 
                failures[i-1] := [];
            end if;
            Append(~failures[i-1], key);

        end for;
    end if;

    prog +:= 1;
    if prog mod 20000 eq 0 then 
        printf "Progress: %o\n", prog;
        for i in Keys(failures) do 
            printf "%o : %o \n", i, #failures[i];
        end for;
        printf "\n";
    end if;
end for;

for i in Keys(failures) do 
    printf "%o : %o \n", i, #failures[i];
end for;
printf "\n";

for i in Keys(failures) do 
    fname := Sprintf("Dataset/zeta_functions/zeta_fails_%o.txt", i);

    str := "";
    failed_keys := failures[i];
    for key in failed_keys do 
        str cat:= Sprint(key) cat "\n";
    end for;
    PrintFile(fname, str);
end for;

F := Open("Dataset/zeta_functions/zeta_fails.m", "w");
WriteObject(F, failures);
delete F;