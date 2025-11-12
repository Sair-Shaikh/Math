


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

ContributionAtFiber := function(q, S, fibA, fibB, fibC)
    if fibA eq 0 then 
        if fibB eq 0 then 
            if fibC eq 0 then 
                return (q+1);
            else 
                return 1;
            end if;
        else 
            if fibC eq 0 then 
                return 2;
            else 
                return 2;
            end if;
        end if;
    else 
        L := fibB/fibA;
        M := fibC/fibA;
    
        if L eq 0 then
            if IsSquare(M) then
                return 1;
            end if;
        elif M/L^2 in S then
            return 2;
        else 
            return 0;
        end if;
    end if;
end function;

PointCount := function(q, A, B, C, A2, B2, C2)
    K := GF(q);
    Embed(k,K);
	S := {x^2 + x : x in K};
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
        
		for z in K do
            //if z eq 0 then continue; end if;
            fibA := Evaluate(NewPolyA, z);
            fibB := Evaluate(NewPolyB, z);
            fibC := Evaluate(NewPolyC, z);
            s1 +:= ContributionAtFiber(q, S, fibA, fibB, fibC)*Osize;
            
		end for;
	end for;

    // Now look at a = 0
    s2 := 0;     
    NewPolyA := Evaluate(A2, [0, 1, r1.1]);
    NewPolyB := Evaluate(B2, [0, 1, r1.1]);
    NewPolyC := Evaluate(C2, [0, 1, r1.1]);

    for i in [1..#FrobeniusReps] do
    	z := FrobeniusReps[i];
		Osize := OrbitSizes[i];

        //if z eq 0 then continue; end if;
        
        fibA := Evaluate(NewPolyA, z);
        fibB := Evaluate(NewPolyB, z);
        fibC := Evaluate(NewPolyC, z);

        s2  +:= ContributionAtFiber(q, S, fibA, fibB, fibC)*Osize;
    end for;

    fibA := Evaluate(A2, [0, 0, 1]);
    fibB := Evaluate(B2, [0, 0, 1]);
    fibC := Evaluate(C2, [0, 0, 1]);

    s3 := 0;
    s3  +:= ContributionAtFiber(q, S, fibA, fibB, fibC);

    printf "s1 = %o, s2 = %o, s3 = %o \n", s1, s2, s3;

    return (s1 + s2 + s3);
end function;



/*
    For each smooth surface X, we want to count points over F_q for q = 2, 4, ..., 2^11. 
    This cell handles the case where X contains a line and has f2 of the form: x1^2 + x2x5 + x3x4. 
*/

F := Open("orbits_with_lines_contain.m", "r");
orbits := ReadObject(F);
delete F;
"Smooth Surfaces Containing a Line: ", #orbits;

F2 := GF(2);
G := GL(5, F2);

R<[x]> := PolynomialRing(F2, 5);
P4<[x]> := ProjectiveSpace(R);

Fiber<a,b,c> := PolynomialRing(F2, 3);
R3<u,v,w>    := PolynomialRing(Fiber, 3);
R2<s, t>     := PolynomialRing(Fiber, 2);

Px<k, l, m> := ProjectiveSpace(F2, 2);

V2, Bit2 := GModule(G, R, 2);
V3, Bit3 := GModule(G, R, 3);

count := 0;

// {orbits[key]["f2int"] : key in Keys(orbits)  };
// {key : key in Keys(orbits) | orbits[key]["f2int"] eq 144 };

for key in Keys(orbits) do

    orbit := orbits[key];

    // Convert from int to polynomials
    f2_int := orbit["f2int"];
    f3_int := orbit["f3int"];
    f2, f3 := convertToPolys(f2_int, f3_int, R, G, Bit2, Bit3);  // Now f2, f3 are in R

    if f2_int eq 144 then continue; end if;
    
    // Pullback f2 to P2[ u : v : w ] with coefficients in F2[a, b, c] -- affine patch not including L
    // other line: q*u + lv*v + lw*w = 0.
    subs_map := hom<R -> R3 | [ u*a, u*b, u*c, v, w]>;
    other_line := subs_map(f2) div u;    
    conic := subs_map(f3) div u;
    
    // Two parametrizations corresponding to q = 0 and q neq 0
    q := MonomialCoefficient(other_line, u);
    lv := MonomialCoefficient(other_line, v);
    lw := MonomialCoefficient(other_line, w);
    
    subs_map_2 := hom<R3 -> R2 | [ lv*s + lw*t, q*s, q*t]>;
    subs_map_3 := hom<R3 -> R2 | [ t, lw*s, lv*s]>;
    other_line := subs_map_2(other_line);
    quadratic := subs_map_2(conic);  // bivariate
    quadratic2 := subs_map_3(conic); // bivariate

    assert other_line eq 0;

    A  := MonomialCoefficient(quadratic, s^2);
    B  := MonomialCoefficient(quadratic, s*t);
    C  := MonomialCoefficient(quadratic, t^2);
    A2 := MonomialCoefficient(quadratic2, s^2);
    B2 := MonomialCoefficient(quadratic2, s*t);
    C2 := MonomialCoefficient(quadratic2, t^2);

    X := Scheme(P4, [f2, f3]);
    time for i in [1..6] do 
        q := 2^i;
        num_pts := PointCount(q, A, B, C, A2, B2, C2);
        printf "q = %o, points = %o \n", q, num_pts;
        
        if q le 64 then 
            Fq := GF(q); 
            XFq := BaseChange(X, Fq);
            n := #Points(XFq);
            printf "Check: q=%o, actual=%o, calc=%o, diff=%o \n", q, n, num_pts, num_pts-n;
        end if;
        
    end for;

    count +:= 1;
    if count ge 10 then 
        break; 
    end if;

end for;

/* 
 656:  c^2*u + b*v + a*w
 1072: b^2*u + c*v + a*w
 1281: a^2*u + c*v + b*w
*/


/*
    For each smooth surface X, we want to count points over F_q for q = 2, 4, ..., 2^11. 
    This cell handles the case where X contains a line and has f2 of the form: x1x4 + x2x5. 
*/

F := Open("orbits_with_lines_contain.m", "r");
orbits := ReadObject(F);
delete F;
"Smooth Surfaces Containing a Line: ", #orbits;

F2 := GF(2);
G := GL(5, F2);

R<[x]> := PolynomialRing(F2, 5);
P4<[x]> := ProjectiveSpace(R);

Fiber<a,b,c> := PolynomialRing(F2, 3);
R3<u,v,w>    := PolynomialRing(Fiber, 3);
R2<s, t>     := PolynomialRing(Fiber, 2);

Px<k, l, m> := ProjectiveSpace(F2, 2);

V2, Bit2 := GModule(G, R, 2);
V3, Bit3 := GModule(G, R, 3);

count := 0;

// {orbits[key]["f2int"] : key in Keys(orbits)  };
// {key : key in Keys(orbits) | orbits[key]["f2int"] eq 144 };

orbits := new_orbits;
for random key in Keys(orbits) do

    orbit := orbits[key];

    // Convert from int to polynomials
    f2_int := orbit["f2int"];
    f3_int := orbit["f3int"];
    f2, f3 := convertToPolys(f2_int, f3_int, R, G, Bit2, Bit3);  // Now f2, f3 are in R

    if not f2_int eq 144 then continue; end if;
    
    // Pullback f2 to P2[ u : v : w ] with coefficients in F2[a, b, c] -- affine patch not including L
    subs_map := hom<R -> R3 | [ u*a, u*b, u*c, v, w]>;
    other_line := subs_map(f2) div u;
    
    conic := subs_map(f3) div u;

    A := Evaluate(MonomialCoefficient(conic, v^2), [0, 0, 1]);
    B := Evaluate(MonomialCoefficient(conic, v*w), [0, 0, 1]);
    C := Evaluate(MonomialCoefficient(conic, w^2), [0, 0, 1]);
    D := Evaluate(MonomialCoefficient(conic, u*v), [0, 0, 1]);
    E := Evaluate(MonomialCoefficient(conic, u*w), [0, 0, 1]);
    F := Evaluate(MonomialCoefficient(conic, u^2), [0, 0, 1]);
    cc := A*l^2 + B*l*m + C*m^2 + D*k*l + E*k*m + F*k^2;
    X_c := Scheme(Px, [cc]);

    q := MonomialCoefficient(other_line, u);
    lv := MonomialCoefficient(other_line, v);
    lw := MonomialCoefficient(other_line, w);

    subs_map_2 := hom<R3 -> R2 | [ t, lw*s, lv*s]>;
    assert subs_map_2(other_line) eq 0;
    conic := subs_map_2(conic);        // should be bivariate now
    
    A := MonomialCoefficient(conic, s^2);
    B := MonomialCoefficient(conic, s*t);
    C := MonomialCoefficient(conic, t^2);

    X := Scheme(P4, [f2, f3]);
    for i in [1..6] do 
        q := 2^i;
        num_pts := PointCount(q, A, B, C, A, B, C);
        Fq := GF(q); 
        X_cq := BaseChange(X_c, Fq);
        pts := #Points(X_cq);
        num_pts +:= pts;
        
        //printf "q = %o, points = %o \n", q, num_pts;
        if q le 64 then 
            Fq := GF(q); 
            XFq := BaseChange(X, Fq);
            n := #Points(XFq);
            printf "Check: q=%o, actual=%o, calc=%o, diff=%o \n", q, n, num_pts, num_pts-n;
        end if;
    end for;

    count +:= 1;
    if count ge 50 then break; end if;
end for;

