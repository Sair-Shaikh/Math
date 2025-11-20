import "magma_scripts/Utils.m" : ConvertToPolys, GetCubicCoeffs;


// Gets counts for f2 orbit neq 144 -- i.e. checks x1^2 + x2x5 + x3x3 not x1x4 + x2x4 and compares against built-in for q <= 64
// Can be used to verify against the C++ counts

k := GF(2);
Fiber<a,b,c> := PolynomialRing(k, 3);
R3<u,s,t> := PolynomialRing(Fiber, 3);

CalculateFrobeniusOrbit := function(K)
	orbits := [];
	for a in K do
		orbit := [a];
		for r in [1..#K] do
			if Frobenius(a,r) ne a then
				Append(~orbit,Frobenius(a,r));
			else
				break;
			end if;
		end for;
	Append(~orbits,Seqset(orbit));
	end for;
	orbits := Setseq(Seqset(orbits));
	orbits := [ Setseq(orbits[i]) : i in [1..#orbits] ];
	FrobeniusReps := [ orbits[i][1] : i in [1..#orbits]];
	OrbitSizes := [ #orbits[i] : i in [1..#orbits]];
    
	return FrobeniusReps, OrbitSizes;
end function;

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

ContributionAtFiber := function(K, A, B, C, D)

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
                return #{x : x in K | x^3+P*x+Q eq 0};        
            end if;
        end if;
    end if;
end function;

PointCount := function(q, A, B, C, D, A2, B2, C2, D2)
    K := GF(q);
    Embed(k,K);
    FrobeniusReps, OrbitSizes := CalculateFrobeniusOrbit(K);
    new_fiber<a, b, c> := ChangeRing(Fiber,K);
    
    r1 := PolynomialRing(K);
	
    // first look at patch where a = 1
	s1 := 0;
	for i in [1..#FrobeniusReps] do
		y := FrobeniusReps[i];
		Osize := OrbitSizes[i];

        NewPolyA := Evaluate(A, [1, y, r1.1]);
        NewPolyB := Evaluate(B, [1, y, r1.1]);
        NewPolyC := Evaluate(C, [1, y, r1.1]);
        NewPolyD := Evaluate(D, [1, y, r1.1]);

		for z in K do
            //if z eq 0 then continue; end if;
            fibA := Evaluate(NewPolyA, z);
            fibB := Evaluate(NewPolyB, z);
            fibC := Evaluate(NewPolyC, z);
            fibD := Evaluate(NewPolyD, z);

            s1 +:= ContributionAtFiber(K, fibA, fibB, fibC, fibD)*Osize;
		end for;
	end for;

    // Now look at a = 0
    s2 := 0;     
    NewPolyA := Evaluate(new_fiber!A2, [0, 1, r1.1]);
    NewPolyB := Evaluate(new_fiber!B2, [0, 1, r1.1]);
    NewPolyC := Evaluate(new_fiber!C2, [0, 1, r1.1]);
    NewPolyD := Evaluate(new_fiber!D2, [0, 1, r1.1]);

    for i in [1..#FrobeniusReps] do
    	z := FrobeniusReps[i];
		Osize := OrbitSizes[i];
        
        fibA := Evaluate(NewPolyA, z);
        fibB := Evaluate(NewPolyB, z);
        fibC := Evaluate(NewPolyC, z);
        fibD := Evaluate(NewPolyD, z);

        s2  +:= ContributionAtFiber(K, fibA, fibB, fibC, fibD)*Osize;
    end for;

    fibA := Evaluate(new_fiber!A2, [0, 0, 1]);
    fibB := Evaluate(new_fiber!B2, [0, 0, 1]);
    fibC := Evaluate(new_fiber!C2, [0, 0, 1]);
    fibD := Evaluate(new_fiber!D2, [0, 0, 1]);

    s3 := 0;
    s3  +:= ContributionAtFiber(K, fibA, fibB, fibC, fibD);

    // printf "s3 = %o, s3+s2 = %o, s1+s2+s3 = %o \n", s3, s2+s3, s1+s2+s3;

    return (s1 + s2 + s3);
end function;

F := Open("Dataset/orbits_lines_secrtpt.m", "r");
orbits := ReadObject(F);
delete F;

"Smooth Surfaces Under Consideration: ", #orbits;
// Set up
F2 := GF(2);
G := GL(5, F2);

R<[x]> := PolynomialRing(F2, 5);
P4<[x]> := ProjectiveSpace(R);

Fiber<a,b,c> := PolynomialRing(F2, 3);
R3<u,v,w>    := PolynomialRing(Fiber, 3);
R2<s, t>     := PolynomialRing(Fiber, 2);

V2, Bit2 := GModule(G, R, 2);
V3, Bit3 := GModule(G, R, 3);

count := 0;
// Main Loop
extensions_range := [1..10];
for key in Keys(orbits) do
    orbit := orbits["5-852295"];
    // orbit := orbits[key];

    // Convert from int to polynomials
    f2_int := orbit["f2int"];
    f3_int := orbit["f3int"];
    f2, f3 := ConvertToPolys(f2_int, f3_int, R, G, Bit2, Bit3);  // Now f2, f3 are in R

    // Pullback f2 to P2[ u : v : w ] with coefficients in F2[a, b, c]
    subs_map := hom<R -> R3 | [ u*a, u*b, u*c, v, w]>;
    other_line := subs_map(f2) div u;    
    elliptic := subs_map(f3);
    
    qq  := MonomialCoefficient(other_line, u);
    lv := MonomialCoefficient(other_line, v);
    lw := MonomialCoefficient(other_line, w);

    if lw eq 0 then // Case 1: lw = 0
        line_para1 := hom<R3 -> R2 | [ lv*s, qq*s,  t]>;
        cubic1  := line_para1(elliptic);  // Deg 3 in P1
        cubic2 := cubic1;

    elif qq eq 0 then // Case 2: qq = 0, lw neq 0
        assert lv eq b and lw eq a;  // Only one orbit of f2 matches this
        line_para1 := hom<R3 -> R2 | [ t, lw*s, lv*s]>;
        cubic1  := line_para1(elliptic);  // Deg 3 in P1
        cubic2  := cubic1;

    else // Case 3: qq neq  0, lw neq 0
        assert lw eq a;
        line_para1 := hom<R3 -> R2 | [ lw*s, lw*t, qq*s+lv*t]>;
        line_para2 := hom<R3 -> R2 | [ lv*s, qq*s, t]>;         
        cubic1  := line_para1(elliptic);   // Deg 3 in P1
        cubic2  := line_para2(elliptic);   // Deg 3 in P1

    end if;

    A, B, C, D := GetCubicCoeffs(cubic1);
    A2, B2, C2, D2 := GetCubicCoeffs(cubic2);

    X := Scheme(P4, [f2, f3]);
    // Now, calculate the correction terms: 
    rt_cubic := Evaluate(f3, [0, 0, 0, s, t]);
    Art, Brt, Crt, Drt := GetCubicCoeffs(rt_cubic);
    for i in extensions_range do 
        q := 2^i;
        Fq := GF(q);
        num_pts := PointCount(q, A, B, C, D, A2, B2, C2, D2);
        
        // Apply corrections
        corr := 0;

        // First, find the number of rational points of X on the line -- should be at least 1
        rt_pts := ContributionAtFiber(Fq, Art, Brt, Crt, Drt);
        assert rt_pts ge 1;
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

        num_pts +:= corr;

        if q le 64 then 
            Fq := GF(q); 
            XFq := BaseChange(X, Fq);
            n := #Points(XFq);
            printf "Check: q=%o, actual=%o, calc=%o, diff=%o \n", q, n, num_pts, num_pts-n;
        else 
            printf "q = %o, points = %o \n", q, num_pts;        
        end if;
    end for;

    count +:= 1;
    if count mod 100 eq 0 then 
        printf "%o/%o\n", count, #orbits; 
    end if;
    break;

end for;