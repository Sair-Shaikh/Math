function ConvertToPolys(f, g, R, G, Bit2, Bit3)
    f := Intseq(f, 2);
    g := Intseq(g, 2);
    f cat:= [0 : i in [1..(15-#f)]];
    f := f @@ Bit2;
    g cat:= [0 : i in [1..(35-#g)]];
    g := g @@ Bit3;
    return f, g;
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


PolyToCpp := function(h)
    if h eq 0 then return "0"; end if;
    if h eq 1 then return "1"; end if;

    ret := "";
    for m in Monomials(h) do
        j := [Degree(m, Parent(h).i) : i in [1..Rank(Parent(h))] ];
        str := "";
        for k in [0 .. #j-1] do
            for l in [1..j[k+1]] do
                if str eq "" then
                    str := Sprintf("y_%o", k);
                else
                    // str := Sprintf("ff2k_mult(%o, y_%o)", str, k);
                    str := Sprintf("mult[%o][y_%o]", str, k);
                end if;
            end for;
        end for;
        
        if ret eq "" then ret := str;
        else ret := ret cat " ^ \\\n   " cat str;
        end if;
    end for;
    return ret;
end function;

CppHeaderTextQuad := function(A, B, C, A2, B2, C2) 
    str := "";
    str cat:= "#define ABC \\" cat "\n";
    str cat:= "    A = " cat PolyToCpp(A) cat ", \\" cat "\n";
    str cat:= "    B = " cat PolyToCpp(B) cat ", \\" cat "\n";
    str cat:= "    C = " cat PolyToCpp(C) cat "\n";

    str cat:= "\n\n";

    str cat:= "#define ABC2 \\" cat "\n";
    str cat:= "    A = " cat PolyToCpp(A2) cat ", \\" cat "\n";
    str cat:= "    B = " cat PolyToCpp(B2) cat ", \\" cat "\n";
    str cat:= "    C = " cat PolyToCpp(C2) cat "\n";

    return str;
end function;


CppHeaderTextCubic := function(A, B, C, D, A2, B2, C2, D2) 
    str := "";
    str cat:= "#define ABCD \\" cat "\n";
    str cat:= "    A = " cat PolyToCpp(A) cat ", \\" cat "\n";
    str cat:= "    B = " cat PolyToCpp(B) cat ", \\" cat "\n";
    str cat:= "    C = " cat PolyToCpp(C) cat ", \\" cat "\n";
    str cat:= "    D = " cat PolyToCpp(D) cat "\n";
    str cat:= "\n\n";

    str cat:= "#define ABCD2 \\" cat "\n";
    str cat:= "    A = " cat PolyToCpp(A2) cat ", \\" cat "\n";
    str cat:= "    B = " cat PolyToCpp(B2) cat ", \\" cat "\n";
    str cat:= "    C = " cat PolyToCpp(C2) cat ", \\" cat "\n";
    str cat:= "    D = " cat PolyToCpp(D2) cat "\n";
    return str;
end function;