import "magma_scripts/Utils.m" : ConvertToPolys, CppHeaderTextQuad;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This script 
// (1) Takes in a dataset containing all K3 surfaces that contain a line x1 = x2 = x3 = 0.
// (2) Parametrizes the quadratic that encodes the roots of the polynomial by projecting from the line
// (3) Writes these as strings in the C++ notation for these functions
// (4) It also computes the deviation the C++ will have from the correct count. This comes from:
//    (4.1) Over counting:  counting points on singular fibers where parametrization is not valid
//    (4.2) Under counting: not counting points on actual singular fibers 
// (5) The corrections are stored with each orbits data. 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

N := 13;
F := Open("Dataset/orbits_lines_contain.m", "r");
orbits := ReadObject(F);
delete F;
"Smooth Surfaces Containing a Line: ", #orbits;

// Clear the coefficients file
fname := "Dataset/cpp_coeffs/quad_coeffs_table_13.txt"; // Writing those that need counting up to 12
F := Open(fname, "w");
delete F; 

// Omit everything that is not solved by coutning up to 13
F := Open("Dataset/zeta_functions/zeta_fails.m", "r");
selected_orbits := ReadObject(F);
delete F;
selected_orbits := selected_orbits[N];

// Setup
F2 := GF(2);
G := GL(5, F2);

R<[x]>       := PolynomialRing(F2, 5);
Fiber<a,b,c> := PolynomialRing(F2, 3);
R3<u,v,w>    := PolynomialRing(Fiber, 3);
R2<s, t>     := PolynomialRing(Fiber, 2);

V2, Bit2 := GModule(G, R, 2);
V3, Bit3 := GModule(G, R, 3);

SetColumns(0);
count := 0;
for key in Keys(orbits) do
    if not key in selected_orbits then continue; end if;

    orbit := orbits[key];

    // Convert from int to polynomials in R<x>
    f2_int := orbit["f2int"];
    f3_int := orbit["f3int"];
    f2, f3 := ConvertToPolys(f2_int, f3_int, R, G, Bit2, Bit3);  // Now f2, f3 are in R
    
    // Pullback f2 to P2[ u : v : w ] with coefficients in F2[a, b, c] and divide by u to work with the blow-up (line is contained with multiplicity 1)
    subs_map := hom<R -> R3 | [ u*a, u*b, u*c, v, w]>;
    other_line := subs_map(f2) div u;    
    conic := subs_map(f3) div u;
    
    // Line N := q * u + lv * v + lw * w = 0. 
    // Kernel contains [lv, q, 0], [lw, 0, q], [0, lw, lv]
    // [lv, q + lw, lv], [lv+lw, q, q], [lw, lw, lv+q] 
    q := MonomialCoefficient(other_line, u);
    lv := MonomialCoefficient(other_line, v);
    lw := MonomialCoefficient(other_line, w);

    if q eq 0 then 
        // This parametrization is valid for all fibers
        assert lv eq b and lw eq a;  // Assumed
        line_para1 := hom<R3 -> R2 | [ t, lw*s, lv*s]>;
        quadratic  := line_para1(conic);  // Deg 2 in P1
        quadratic2 := quadratic;

    else // q neq  0
        // We need two parametrizations for q = 0 and q neq 0. 
        // We assume surfaces are such that if q \neq 0, then q = a^2
        assert q eq a^2;
        line_para1 := hom<R3 -> R2 | [ lv*s + lw*t, q*s, q*t]>; // Span([lv, 0, q], [lw, q, 0])) -- when q neq 0, this is valid
        line_para2 := hom<R3 -> R2 | [ t, lw*s, lv*s]>;         // Span([1, 0, 0],  [0, lw, lv]) -- when q = 0, this is valid

        quadratic  := line_para1(conic);  // Deg 2 in P1
        quadratic2 := line_para2(conic);  // Deg 2 in P1

    end if;


    // Now, calculate the correction terms: 

    corrections := [];
    if q eq 0 then 
        // In this case, the parametrization vanishes over the final fiber.
        // We plug in [t, 0, 0].
        // The quadratic is then t^2, if the conic had a u^2 term, else 0.
        // We overcount by 1 or by q+1. 

        // We need to count points on the conic to correct the counts
        A := Evaluate(MonomialCoefficient(conic, v^2), [0, 0, 1]);
        B := Evaluate(MonomialCoefficient(conic, v*w), [0, 0, 1]);
        C := Evaluate(MonomialCoefficient(conic, w^2), [0, 0, 1]);
        D := Evaluate(MonomialCoefficient(conic, u*v), [0, 0, 1]);
        E := Evaluate(MonomialCoefficient(conic, u*w), [0, 0, 1]);
        F := Evaluate(MonomialCoefficient(conic, u^2), [0, 0, 1]);

        HasU2Term := Integers()!F; // 1 if yes, 0 if not

        cc := A*b^2 + B*b*c + C*c^2 + D*a*b + E*a*c + F*a^2; // reusing a, b, c -- doesnt matter
        disc := A*E^2 + B^2*F + C*D^2 -B*D*E; 
        
        for i in [1..N] do 
            q := 2^i;
            over  := (1-HasU2Term)*q+1;
            if not disc eq 0 then 
                // Conic is smooth -- has q+1 points
                under := (q+1);
            else
                // TODO: Conic is not smooth -- could do Arf Invariants
                K := GF(q);
                under := #Points(Curve(Proj(Fiber), cc), K);
            end if;
            Append(~corrections, under-over);
        end for;
    else 
        corrections := [ 0 : i in [1..N]];
    end if;


    // Write coefficients and corrections to text file -- we will import this into Python to do C++ orchestration
    A  := MonomialCoefficient(quadratic, s^2);
    B  := MonomialCoefficient(quadratic, s*t);
    C  := MonomialCoefficient(quadratic, t^2);
    A2 := MonomialCoefficient(quadratic2, s^2);
    B2 := MonomialCoefficient(quadratic2, s*t);
    C2 := MonomialCoefficient(quadratic2, t^2);

    str := "Key: " cat key cat "\n";
    str cat:= CppHeaderTextQuad(A, B, C, A2, B2, C2);
    str cat:= "\n";
    str cat:= "Corrections: " cat Sprint(corrections);
    PrintFile(fname, str);

    count +:= 1;
    if count mod 2000 eq 0 then 
        printf "Progress: %o/%o\n", count, #orbits;
    end if;

end for;




