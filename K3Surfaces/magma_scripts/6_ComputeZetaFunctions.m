
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

    Using the functional equation, we can determine #X(F_q^n) for n = 12, ..., 22 using #X(F_q^n) for n = 1, ,,, 11.
    (Details)

    Using Newton's formulas


*/


HalfWeil := function(point_counts) 

    assert #point_counts eq 11;

    traces := [(point_counts[i] - 1 - 2^(2*i)) : i in [1..11]];
    // c := [Rationals() ! -2^31 : i in [1..11]]; // Initialize
    // c[1] := -traces[1];
    // for k in [2..11] do
    //     c[k] := (traces[k] + &+[c[i]*traces[k-i] : i in [1..k-1]])/(-k);
    // end for;
    // [1] cat c;
    // FrobeniusTracesToWeilPolynomials(traces, 2, 2, 22)[1];
    return FrobeniusTracesToWeilPolynomials(traces, 2, 2, 22);
end function;


F := Open("Dataset/orbits_point_counts.m", "r");
orbits := ReadObject(F);
delete F;

counts := AssociativeArray();
counts[0] := 0;
counts[1] := 0;
counts[2] := 0;

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
    f2, f3 := ConvertToPolys(f2_int, f3_int, R, G, Bit2, Bit3);  // Now f2, f3 are in R

    pt_count := orbit["pt_count"];
    X := Scheme(ProjectiveSpace(R), [f2, f3]);

    for i in [1..4] do 
        q := 2^i;
        Fq := GF(q); 
        XFq := BaseChange(X, Fq);
        n := #Points(XFq);

        if not n eq pt_count[i] then 
            key;
            break;
        end if;

    end for;
    
    char_polys := HalfWeil(pt_count);

    if #char_polys eq 2 then 
        // Find lowest index above 11 for which coeff is 0
        coeffs := Coefficients(char_polys[1]);
        assert coeffs[12] eq 0;

        for i in [1..11] do 
            if not coeffs[12+i] eq 0 then 
                if not IsDefined(failures, i) then 
                    failures[i] := 0;
                end if;
                failures[i] +:= 1; 
                break;
            end if;
        end for;

    end if;

    prog +:= 1;
    if prog mod 1000 eq 0 then 
        printf "%o\n", prog;
        for i in Keys(failures) do 
            printf "%o : %o \n", 11+i, failures[i];
        end for;
        printf "\n";
    end if;


end for;
