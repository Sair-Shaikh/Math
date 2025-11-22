function ConvertToPolys(f, g, R, G, Bit2, Bit3)
    f := Intseq(f, 2);
    g := Intseq(g, 2);
    f cat:= [0 : i in [1..(15-#f)]];
    f := f @@ Bit2;
    g cat:= [0 : i in [1..(35-#g)]];
    g := g @@ Bit3;
    return f, g;
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
                return #{e : e in K | e^3 + K!P*e + K!Q eq 0};        
            end if;
        end if;
    end if;
end function;


GetCubicCoeffs := function(f)
    R := Parent(f);
    t := R.1;
    s := R.2;

    A := MonomialCoefficient(f, t^3);
    B := MonomialCoefficient(f, t^2*s);
    C := MonomialCoefficient(f, t*s^2);
    D := MonomialCoefficient(f, s^3);

    // Depressed cubic coefficients t^3 + Pts^2 + Qs^3 -- substituting t -> t+B
    if A eq 1 then 
        P := B^2+C;
        Q := B*C+D;
        return A, 0, P, Q;
    elif D eq 1 then 
        P := C^2+B;
        Q := B*C+A;
        return D, 0, P, Q;
    else 
        return A, B, C, D;
    end if;

end function;

PolyToCpp := function(h, typ)
    if h eq 0 then return "0"; end if;
    if h eq 1 then return "1"; end if;

    ret := "";
    for m in Monomials(h) do
        j := [Degree(m, Parent(h).i) : i in [1..Rank(Parent(h))] ];
        str := "";
        for k in [0 .. #j-1] do
            for l in [1..j[k+1]] do
                if str eq "" then
                    str := Sprintf("y_%o", k-1);
                else
                    if typ eq 1 then 
                        str := Sprintf("mult[%o][y_%o]", str, k);
                    elif typ eq 2 then 
                        str := Sprintf("ff2k_mult(%o, y_%o)", str, k-1);
                    else
                        assert false;
                    end if;
                end if;
            end for;
        end for;
        
        if ret eq "" then ret := str;
        else ret := ret cat " ^ \\\n   " cat str;
        end if;
    end for;
    return ret;
end function;


PolyToBitVector := function(f)
    // if f eq 0 then return "\"0\""; end if;
    // if f eq 1 then return "\"1\""; end if;
    // return "\"" cat Sprint(Seqint(bit_array, 2)) cat "\""; // Returns an integer whos bit value is the required
    
    if f eq 0 then
        set_mons := [0] cat [0 : i in [1..29]];    
    elif f eq 1 then 
        set_mons := [1] cat [0 : i in [1..29]];
    else     
        assert Degree(f) le 6;
        R<a,b,c> := PolynomialRing(GF(2), 3, "lex");
        f_h := Generators(Homogenization(ideal< R | [f]>, true, "lex"))[1];
        Rh := Parent(f_h);
        f_h := Rh.1^(6-Degree(f_h)) * f_h; // Rh.1 is the homogenizing variable
        bit_array := [ Integers()!MonomialCoefficient(f_h, mon)  : mon in MonomialsOfDegree(Rh, 6)];
        set_mons := [ i : i in [1..#bit_array] | bit_array[i] eq 1];
        set_mons cat:= [0 : i in [1..(30-#set_mons)]];
    end if;

    str_rep := Sprint(set_mons); 
    str_rep := "{ " cat str_rep[2..(#str_rep-1)] cat "}";

    return str_rep;
end function;



CppHeaderTextQuad := function(A, B, C, A2, B2, C2) 
    str := "";
    str cat:= "#define ABC \\" cat "\n";
    str cat:= "    A = " cat PolyToCpp(A, 1) cat ", \\" cat "\n";
    str cat:= "    B = " cat PolyToCpp(B, 1) cat ", \\" cat "\n";
    str cat:= "    C = " cat PolyToCpp(C, 1) cat "\n";

    str cat:= "\n\n";

    str cat:= "#define ABC2 \\" cat "\n";
    str cat:= "    A = " cat PolyToCpp(A2, 1) cat ", \\" cat "\n";
    str cat:= "    B = " cat PolyToCpp(B2, 1) cat ", \\" cat "\n";
    str cat:= "    C = " cat PolyToCpp(C2, 1) cat "\n";

    return str;
end function;


CppHeaderTextCubic := function(A, B, C, D, A2, B2, C2, D2) 
    str := "";
    str cat:= "#define ABC \\" cat "\n";
    str cat:= "    A = " cat PolyToCpp(A, 1) cat ", \\" cat "\n";
    str cat:= "    B = " cat PolyToCpp(B, 1) cat ", \\" cat "\n";
    str cat:= "    C = " cat PolyToCpp(C, 1) cat ", \\" cat "\n";
    str cat:= "    D = " cat PolyToCpp(D, 1) cat "\n";
    str cat:= "\n\n";

    str cat:= "#define ABC2 \\" cat "\n";
    str cat:= "    A = " cat PolyToCpp(A2, 1) cat ", \\" cat "\n";
    str cat:= "    B = " cat PolyToCpp(B2, 1) cat ", \\" cat "\n";
    str cat:= "    C = " cat PolyToCpp(C2, 1) cat ", \\" cat "\n";
    str cat:= "    D = " cat PolyToCpp(D2, 1) cat "\n";
    return str;
end function;


CppHeaderTextQuadBigQ := function(A, B, C, A2, B2, C2) 
    str := "";
    str cat:= "#define ABC \\" cat "\n";
    str cat:= "    A = " cat PolyToCpp(A, 2) cat ", \\" cat "\n";
    str cat:= "    B = " cat PolyToCpp(B, 2) cat ", \\" cat "\n";
    str cat:= "    C = " cat PolyToCpp(C, 2) cat "\n";

    str cat:= "\n\n";

    str cat:= "#define ABC2 \\" cat "\n";
    str cat:= "    A = " cat PolyToCpp(A2, 2) cat ", \\" cat "\n";
    str cat:= "    B = " cat PolyToCpp(B2, 2) cat ", \\" cat "\n";
    str cat:= "    C = " cat PolyToCpp(C2, 2) cat "\n";

    return str;
end function;

CppHeaderTextCubicBigQ := function(A, B, C, D, A2, B2, C2, D2) 
    str := "";
    str cat:= "#define ABC \\" cat "\n";
    str cat:= "    A = " cat PolyToCpp(A, 2) cat ", \\" cat "\n";
    str cat:= "    B = " cat PolyToCpp(B, 2) cat ", \\" cat "\n";
    str cat:= "    C = " cat PolyToCpp(C, 2) cat ", \\" cat "\n";
    str cat:= "    D = " cat PolyToCpp(D, 2) cat "\n";
    str cat:= "\n\n";

    str cat:= "#define ABC2 \\" cat "\n";
    str cat:= "    A = " cat PolyToCpp(A2, 2) cat ", \\" cat "\n";
    str cat:= "    B = " cat PolyToCpp(B2, 2) cat ", \\" cat "\n";
    str cat:= "    C = " cat PolyToCpp(C2, 2) cat ", \\" cat "\n";
    str cat:= "    D = " cat PolyToCpp(D2, 2) cat "\n";
    return str;
end function;

// CppHeaderTextCubicFast := function(A, B, C, D, A2, B2, C2, D2) 
//     str := "";
//     str cat:= "#define ABC \\" cat "\n";
//     str cat:= "    A = to_uint128(" cat PolyToBitVector(A) cat "), \\" cat "\n";
//     str cat:= "    B = to_uint128(" cat PolyToBitVector(B) cat "), \\" cat "\n";
//     str cat:= "    C = to_uint128(" cat PolyToBitVector(C) cat "), \\" cat "\n";
//     str cat:= "    D = to_uint128(" cat PolyToBitVector(D) cat ")\n";
//     str cat:= "\n\n";

//     str cat:= "#define ABC2 \\" cat "\n";
//     str cat:= "    A = to_uint128(" cat PolyToBitVector(A2) cat "), \\" cat "\n";
//     str cat:= "    B = to_uint128(" cat PolyToBitVector(B2) cat "), \\" cat "\n";
//     str cat:= "    C = to_uint128(" cat PolyToBitVector(C2) cat "), \\" cat "\n";
//     str cat:= "    D = to_uint128(" cat PolyToBitVector(D2) cat ")\n";
//     return str;
// end function;


// CppHeaderTextQuadFast := function(A, B, C, A2, B2, C2) 
//     str := "";
//     str cat:= "#define ABC \\" cat "\n";
//     str cat:= "    A = to_uint128(" cat PolyToBitVector(A) cat "), \\" cat "\n";
//     str cat:= "    B = to_uint128(" cat PolyToBitVector(B) cat "), \\" cat "\n";
//     str cat:= "    C = to_uint128(" cat PolyToBitVector(C) cat ")\n";
//     str cat:= "\n";

//     str cat:= "#define ABC2 \\" cat "\n";
//     str cat:= "    A = to_uint128(" cat PolyToBitVector(A2) cat "), \\" cat "\n";
//     str cat:= "    B = to_uint128(" cat PolyToBitVector(B2) cat "), \\" cat "\n";
//     str cat:= "    C = to_uint128(" cat PolyToBitVector(C2) cat ")\n";
//     return str;
// end function;

CppHeaderTextCubicFast := function(A, B, C, D, A2, B2, C2, D2) 
    str := "";
    str cat:= "#define ABC \\" cat "\n";
    str cat:= "    A = " cat PolyToBitVector(A) cat ", \\" cat "\n";
    str cat:= "    B = " cat PolyToBitVector(B) cat ", \\" cat "\n";
    str cat:= "    C = " cat PolyToBitVector(C) cat ", \\" cat "\n";
    str cat:= "    D = " cat PolyToBitVector(D) cat "\n";
    str cat:= "\n\n";

    str cat:= "#define ABC2 \\" cat "\n";
    str cat:= "    A = " cat PolyToBitVector(A2) cat ", \\" cat "\n";
    str cat:= "    B = " cat PolyToBitVector(B2) cat ", \\" cat "\n";
    str cat:= "    C = " cat PolyToBitVector(C2) cat ", \\" cat "\n";
    str cat:= "    D = " cat PolyToBitVector(D2) cat "\n";
    return str;
end function;


CppHeaderTextQuadFast := function(A, B, C, A2, B2, C2) 
    str := "";
    str cat:= "#define ABC \\" cat "\n";
    str cat:= "    A = " cat PolyToBitVector(A) cat ", \\" cat "\n";
    str cat:= "    B = " cat PolyToBitVector(B) cat ", \\" cat "\n";
    str cat:= "    C = " cat PolyToBitVector(C) cat "\n";
    str cat:= "\n";

    str cat:= "#define ABC2 \\" cat "\n";
    str cat:= "    A = " cat PolyToBitVector(A2) cat ", \\" cat "\n";
    str cat:= "    B = " cat PolyToBitVector(B2) cat ", \\" cat "\n";
    str cat:= "    C = " cat PolyToBitVector(C2) cat "\n";
    return str;
end function;