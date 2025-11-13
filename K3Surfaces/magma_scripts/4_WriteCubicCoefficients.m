import "magma_scripts/Utils.m" : ConvertToPolys, CppHeaderTextCubic, GetCubicCoeffs;

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

NumQuadSols := function(K, A, B, C)
    // This only counts affine solutions 
    if A eq 0 then 
        if B eq 0 then 
            if C eq 0 then 
                return #K;
            else 
                return 0;
            end if;
        else 
            return 1;
        end if;
    else 
        L := B/A;
        M := C/A;
    
        if L eq 0 then
            return 1;
        elif M/L^2 in {x^2 + x : x in K} then
            return 2;
        else 
            return 0;
        end if;
    end if;

end function;

NumCubicSols := function(K, A, B, C, D)

    if A eq 0 then 
        // This is at most a quadratic -- solve that
        // If [ t : 1], then get NumSols of quadratic. 
        // If [ t : 0], then the equation vanishes so count that
        return 1+NumQuadSols(K, B, C, D); 

    else 
        // Make it depressed:
        P := (B/A)^2 + (C/A);
        Q := (B/A)*(C/A) + (D/A);


        if Q eq 0 then 
            if P eq 0 then  // t^3 = 0
                return 1;
            else            // t^3 + Pt = 0 -> t = 0, t = \sqrt{P}
                if IsSquare(P) then
                    return 2;
                else 
                    return 1;
                end if;
            end if;
        else 
            if P eq 0 then // t^3 + Q = 0
                // This part is not correct
                if Q in {e^3 : e in K} then
                    if (#K-1) mod 3 eq 0 then 
                        return 3;
                    else
                        return 1;
                    end if;
                else       
                    return 0;
                end if;
            else          // t^3 + Pt + Q = 0
                // Just default to looping -- no need to Cardano's for such few solves
                P; Q;
                return #{e : e in K | e^3 + K!P*e + K!Q eq 0};        
            end if;
        end if;
    end if;
end function;

F := Open("Dataset/orbits_lines_sec.m", "r");
orbits := ReadObject(F);
delete F;
"Smooth Surfaces With A Secant Line: ", #orbits;

batch := 1;
fname := Sprintf("Dataset/CppCoeffs/container_cube_coeffs_table_1%o.txt", batch);

// Clear the coefficients file
F := Open(fname, "w");
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

new_orbits := AssociativeArray();
count := 0;
time for key in Keys(orbits) do

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
    for i in [1..11] do 
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
    end for;

    // Write coefficients and corrections to text file -- we will import this into Python to do C++ orchestration
    str := "Key: " cat key cat "\n";
    str cat:= CppHeaderTextCubic(A, B, C, D, A2, B2, C2, D2);
    str cat:= "\n";
    str cat:= "Corrections: " cat Sprint(corrections);
    PrintFile(fname, str);

    // Also save these in a magma file
    new_orbits[key] := orbit;
    new_orbits[key]["ptct_corr"] := corrections;

    count +:= 1;
    if count mod 1000 eq 0 then 
        printf "Progress: %o/%o\n", count, #orbits;
    end if;


    if count mod 100000 eq 0 then 
        // Change the coefficients File
        batch +:= 1;
        fname := Sprintf("Dataset/CppCoeffs/container_cube_coeffs_table_%o.txt", batch);

        // Clear the new coefficients file
        F := Open(fname, "w");
        delete F; 

    end if;

end for;


F := Open("Dataset/orbits_count_corrs_secanters.m", "r");
old_orbits := ReadObject(F);
delete F;

#old_orbits;
for k -> v in new_orbits do 
    old_orbits[k] := v;
end for;
#old_orbits;

F := Open("Dataset/orbits_count_corrs_secanters2.m", "w");
WriteObject(F, old_orbits);
delete F;



