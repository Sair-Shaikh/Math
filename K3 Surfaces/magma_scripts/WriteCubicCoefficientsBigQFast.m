import "magma_scripts/Utils.m" : ConvertToPolys, CppHeaderTextCubicFast, GetCubicCoeffs, NumQuadSols, NumCubicSols;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This script 
// (1) Takes in a dataset containing all K3 surfaces that have a secant line x1 = x2 = x3 = 0, contained in V(f2)
// (2) Parametrizes the cubic that encodes the roots of the polynomial by projecting from the line
// (3) Writes these as strings in the C++ notation for these functions
// (4) It also computes the deviation the C++ will have from the correct count. This comes from:
//    (4.1) Over counting:  counting points on singular fibers where parametrization is not valid
//    (4.2) Under counting: not counting points on actual singular fibers
//    (4.3) The blowup has a P1 for each rational point
// (5) The corrections are stored with each orbits data. 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


F := Open("Dataset/orbits_count_corrs_secanters2.m", "r");
orbits := ReadObject(F);
delete F;
"Smooth Surfaces With A Secant Line: ", #orbits;

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
RP1<y,z>     := PolynomialRing(F2, 2);

V2, Bit2 := GModule(G, R, 2);
V3, Bit3 := GModule(G, R, 3);
SetColumns(0);

for N in [14..20] do 
    // Select only those need counting up to N
    selected_orbits := all_selected_orbits[N];
    fname := Sprintf("Dataset/cpp_coeffs/cubic_bigq_%o_fast.txt", N);
    printf "Selected Orbits: %o\n", #selected_orbits;

    if #selected_orbits eq 0 then continue; end if;

    // Clear the coefficients file
    F := Open(fname, "w");
    delete F; 

    count := 0;

    time for key in selected_orbits do
        
        if not IsDefined(orbits, key) then continue; end if; // key might refer to surfaces with a line

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
        
        // Pullback f2 to P2[ u : v : w ] with coefficients in F2[a, b, c] and divide by u to work with the blow-up
        subs_map := hom<R -> R3 | [ u*a, u*b, u*c, v, w]>;
        other_line := subs_map(f2) div u;    

        // Also pullback f3 to (generically) an elliptic curve
        elliptic := subs_map(f3);
        
        // Line N := q * u + lv * v + lw * w = 0
        qq  := MonomialCoefficient(other_line, u);
        lv := MonomialCoefficient(other_line, v);
        lw := MonomialCoefficient(other_line, w);


        if lw eq 0 then // Case 1: lw = 0
            // This parametrization is valid for all fibers -- except the final one
            line_para1 := hom<R3 -> R2 | [ lv*s, qq*s,  t]>;    // Span([lv, qq, 0],  [0, 0, 1]) -- when lw = 0, this is valid
            cubic1  := line_para1(elliptic);  // Deg 3 in P1
            cubic2 := cubic1;

        elif qq eq 0 then // Case 2: qq = 0, lw neq 0
            // This parametrization is valid for all fibers -- except the final one
            assert lv eq b and lw eq a;  // Only one orbit of f2 matches this
            line_para1 := hom<R3 -> R2 | [ t, lw*s, lv*s]>;
            cubic1  := line_para1(elliptic);  // Deg 3 in P1
            cubic2  := cubic1;

        else // Case 3: qq neq  0, lw neq 0
            // By quirks of how we picked orbit representatives, we minimized degree in x[5]
            // Thus, very few of these have a lw -- those that have have been made to have lw = a, making fibration easier
            assert lw eq a;
            line_para1 := hom<R3 -> R2 | [ lw*s, lw*t, qq*s+lv*t]>; // Span([lw, 0, qq], [0, lw, lv])) -- when lw neq 0, this is valid
            line_para2 := hom<R3 -> R2 | [ lv*s, qq*s, t]>;         // Span([0, 0, 1],  [lv, qq, 0])   -- when lw = 0, this is valid
            cubic1  := line_para1(elliptic);   // Deg 3 in P1
            cubic2  := line_para2(elliptic);   // Deg 3 in P1

        end if;

        A, B, C, D := GetCubicCoeffs(cubic1);
        A2, B2, C2, D2 := GetCubicCoeffs(cubic2);
        
        // Now, calculate the correction terms: 
        rt_cubic := Evaluate(f3, [0, 0, 0, s, t]);
        Art, Brt, Crt, Drt := GetCubicCoeffs(rt_cubic);
        corrections := [];

        i := N;
        q := 2^i;
        Fq := GF(q);
        Embed(F2, Fq);
        corr := 0;

        // First, find the number of rational points of X on the line
        rt_pts := NumCubicSols(Fq, Art, Brt, Crt, Drt);
        corr +:= -q*rt_pts;

        // Contribution over singular fibers, i.e. where f2 vanishes -- do this in the dumbest way possible too
        P2Fq       := ProjectiveSpace(Fq, 2);
        R3Fq       := ChangeRing(R3, Fq);
        singular_points := Points(Scheme(P2Fq, [qq, lv, lw]), Fq);
        coeffs, mons := CoefficientsAndMonomials(elliptic);
        for fib in singular_points do 
            new_coeffs := [ Evaluate(co, Eltseq(fib)) : co in coeffs ];
            f3_new     := &+[ new_coeffs[i]*R3Fq!mons[i] : i in [1..#new_coeffs] ];
            corr      +:= #Points(Curve(P2Fq, f3_new))-1;  // -1 because the "incorrect" fibration will count 1 point
        end for;
        Append(~corrections, corr);

        // Write coefficients and corrections to text file -- we will import this into Python to do C++ orchestration
        str := "Key: " cat key cat "\n";
        str cat:= CppHeaderTextCubicFast(A, B, C, D, A2, B2, C2, D2);
        str cat:= "\n";
        str cat:= "Corrections: " cat Sprint(corrections);
        PrintFile(fname, str);

        count +:= 1;
        if count mod 1000 eq 0 then 
            printf "Progress: %o/%o\n", count, #orbits;
        end if;

    end for;
end for;