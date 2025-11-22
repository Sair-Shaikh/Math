import "magma_scripts/Utils.m" : ConvertToPolys, CppHeaderTextQuadFast;

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

F := Open("Dataset/orbits_lines_contain.m", "r");
orbits := ReadObject(F);
delete F;
"Smooth Surfaces Containing a Line: ", #orbits;

// Select only those that need additional counting
F := Open("Dataset/zeta_functions/zeta_fails.m", "r");
all_selected_orbits := ReadObject(F);
delete F;

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


for N in [14..20] do
    // Select only those that need counting up to N
    if not IsDefined(all_selected_orbits, N) then continue; end if;
    selected_orbits := all_selected_orbits[N];
    printf "Selected Orbits: %o\n", #selected_orbits;

    if #selected_orbits eq 0 then continue; end if;

    // Clear the coefficients file
    fname := Sprintf("Dataset/cpp_coeffs/quad_coeffs_%o_fast_test.txt", N); // Writing those that need counting up to 12
    F := Open(fname, "w");
    delete F; 

    count := 0;
    for key in selected_orbits do
        if not IsDefined(orbits, key) then continue; end if; // key might refer to surfaces with a line

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
        // Kernel contains [lv, qq, 0], [lw, 0, qq], [0, lw, lv]
        // [lv, qq + lw, lv], [lv+lw, qq, qq], [lw, lw, lv+qq] 
        qq := MonomialCoefficient(other_line, u);
        lv := MonomialCoefficient(other_line, v);
        lw := MonomialCoefficient(other_line, w);

        if qq eq 0 then 
            // This parametrization is valid for all fibers
            assert lv eq b and lw eq a;  // Assumed
            line_para1 := hom<R3 -> R2 | [ t, lw*s, lv*s]>;
            quadratic  := line_para1(conic);  // Deg 2 in P1
            quadratic2 := quadratic;

        else // qq neq  0
            // We need two parametrizations for qq = 0 and qq neq 0. 
            // We assume surfaces are such that if qq \neq 0, then qq = a^2
            assert qq eq a^2;
            line_para1 := hom<R3 -> R2 | [ lv*s + lw*t, qq*s, qq*t]>; // Span([lv, 0, qq], [lw, qq, 0])) -- when qq neq 0, this is valid
            line_para2 := hom<R3 -> R2 | [ t, lw*s, lv*s]>;         // Span([1, 0, 0],  [0, lw, lv]) -- when qq = 0, this is valid

            quadratic  := line_para1(conic);  // Deg 2 in P1
            quadratic2 := line_para2(conic);  // Deg 2 in P1

        end if;
        // Now, calculate the correction terms: 

        corrections := [];
        if qq eq 0 then 
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
            
            q := 2^N;
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
        else 
            corrections := [ 0 ];
        end if;


        // Write coefficients and corrections to text file -- we will import this into Python to do C++ orchestration
        A  := MonomialCoefficient(quadratic, s^2);
        B  := MonomialCoefficient(quadratic, s*t);
        C  := MonomialCoefficient(quadratic, t^2);
        A2 := MonomialCoefficient(quadratic2, s^2);
        B2 := MonomialCoefficient(quadratic2, s*t);
        C2 := MonomialCoefficient(quadratic2, t^2);

        str := "Key: " cat key cat "\n";
        str cat:= CppHeaderTextQuadFast(A, B, C, A2, B2, C2);
        str cat:= "Corrections: " cat Sprint(corrections) cat "\n";
        PrintFile(fname, str);

        count +:= 1;
        if count mod 2000 eq 0 then 
            printf "Progress: %o/%o\n", count, #orbits;
        end if;

    end for;
end for;