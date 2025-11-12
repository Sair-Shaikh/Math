k  := FiniteField(2);
k2 := VectorSpace(k, 2);
k5 := VectorSpace(k, 5);
R<[x]> := PolynomialRing(k, 5);
P4<[x]> := ProjectiveSpace(R);
G := GL(5, k);
V2, Bit2 := GModule(G, R, 2);
V3, Bit3 := GModule(G, R, 3);

function convertToPolys(f, g, R, G, Bit2, Bit3)
    f := Intseq(f, 2);
    g := Intseq(g, 2);
    f cat:= [0 : i in [1..(15-#f)]];
    f := f @@ Bit2;
    g cat:= [0 : i in [1..(35-#g)]];
    g := g @@ Bit3;
    return f, g;
end function;

function convertLineToStandardForm(f2, f3, rref) 

    plane := Image(rref);
    complement := Complement(k5, plane);

    line_basis := Basis(complement);

    new_basis := line_basis cat [rref[1], rref[2]];
    inv_change_of_basis := Transpose(Matrix(5, 5, new_basis));  // Matrix with columns = this basis
    subs  := [ &+[ inv_change_of_basis[i,j] * x[j] : j in [1..5] ] : i in [1..5] ];

    f2 := Evaluate(f2, subs);
    f3 := Evaluate(f3, subs);  
    
    
    // Sanity checks: 
    /*
    change_of_basis := inv_change_of_basis^-1;
    line   := Transpose(change_of_basis * Transpose(Matrix(2, 5, [rref[1], rref[2]])));
    assert line eq Matrix(2, 5, [ [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]]);
    */

    return f2, f3;

end function;


function LineIsContained(f2, f3, rref)
    R2<u,v> := PolynomialRing(k, 2);
    subs := [rref[1][i]*u + rref[2][i]*v : i in [1..5]];
    f2_subbed := Evaluate(f2, subs);
    f3_subbed := Evaluate(f3, subs);

    return (IsZero(f2_subbed) and IsZero(f3_subbed));
end function;


F := Open("orbit_reps_smooth.m", "r");
K3_smooth := ReadObject(F);
delete F;
"Total number of smooth surfaces: ", #K3_smooth;


/* 
    For each surface, sort into:
        (1) Surfaces containing a line.
        (2) Surfaces with a k-secant line, for k = 2, 3, 4, 5, 6
    Convert the polynomials to make the line be in standard format. 
*/

k  := FiniteField(2);
k2 := VectorSpace(k, 2);
k5 := VectorSpace(k, 5);

// All lines in P4, i.e. points of Grassmanian(2, 5)
echforms := SetToSequence({EchelonForm(M) : M in Hom(k2, k5) | Rank(M) eq 2}); // 155 of these exist
"Total number of lines: ", #echforms;

G := GL(5, k);
R<[x]> := PolynomialRing(k, 5);
R5<a, b, c, s, t> := PolynomialRing(k, 5);
P4<[x]> := ProjectiveSpace(R);
R2<u,v> := PolynomialRing(k, 2);

V2, Bit2 := GModule(G, R, 2);
V3, Bit3 := GModule(G, R, 3);

orbs_with_lines := AssociativeArray();
containers := AssociativeArray();
rat_pters := AssociativeArray();
leftovers := AssociativeArray();


unique_f2s := [];
progress := 0;
for key in Keys(K3_smooth) do
    orbit := K3_smooth[key];
    f2_int := orbit["f2int"];
    f3_int := orbit["f3int"];
    f2, f3 := convertToPolys(f2_int, f3_int, R, G, Bit2, Bit3);

    lines_data := AssociativeArray();
    contains_line := false;
    secant_with_ratpt := false;
    max_secancy := 0;
    for rref in echforms do 

        // Complete classification of lines in each surface
        //   (1) Line is contained
        //   (2) Line is k-secant and has rational points
        //   (3) Line is k-secant and has no rational points

        // Pullback f2 and f3 to the line
        subs := [rref[1][i]*u + rref[2][i]*v : i in [1..5]];
        f2_subbed := Evaluate(f2, subs);
        f3_subbed := Evaluate(f3, subs);

        // Check if the line is contained
        if (IsZero(f2_subbed) and IsZero(f3_subbed)) then
            
            if not IsDefined(lines_data, -1) then lines_data[-1] := 0; end if;
            lines_data[-1] +:= 1;

            if not contains_line then // have not already found a contained line
                selected_line := rref;
                contains_line := true;
            end if;
            
            continue;
        end if;

        // Calculate the secancy of the line
        secancy := Degree(GreatestCommonDivisor(f2_subbed, f3_subbed));

        // Calculate if current line has rational pts over F2
        cur_rt_pt := false;
        pts := [ [0, 1], [1, 0], [1, 1] ]; 

        for pt in pts do 
            if IsZero(Evaluate(f2_subbed, pt)) and IsZero(Evaluate(f3_subbed, pt)) then 
                cur_rt_pt := true; 
                break;
            end if;
        end for; 
        
        if not contains_line and secancy ge max_secancy then 

            if (not secant_with_ratpt or cur_rt_pt) then 
                secant_with_ratpt := cur_rt_pt; 
                selected_line := rref;
            end if;

            max_secancy := secancy;
            
        end if;

        // Add the line data to lines_data
        if not IsDefined(lines_data, secancy) then lines_data[secancy] := 0; end if;
        lines_data[secancy] +:= 1;

    end for;
    
    f2, f3 := convertLineToStandardForm(f2, f3, selected_line);
    orbit["f2int"] := Seqint([Integers() !Bit2(f2)[i] : i in [1..#Basis(V2)]], 2);
    orbit["f3int"] := Seqint([Integers() !Bit3(f3)[i] : i in [1..#Basis(V3)]], 2);
    orbit["lines_data"] := lines_data;
    orbs_with_lines[key] := orbit;

    // Maintain counts and partition into the correct dicts
    if contains_line then 
        containers[key] := orbit;
    elif secant_with_ratpt then 
        rat_pters[key] := orbit;
    else 
        leftovers[key] := orbit;
    end if;
    
    progress +:= 1;
    if progress mod 1000 eq 0 then
        "Progress: ", progress;
        "Contains Line: ", #containers;
        "Sec with RatPt: ", #rat_pters;
        "The Rest: ", #leftovers;
        "";
    end if;
end for;



/*
    Post-Process by permuting x1, x2, x3 to be in standard forms. 
*/

orbits := containers;
new_orbits := AssociativeArray();
count := 0;
for key in Keys(orbits) do


    orbit := orbits[key];
    
    // Convert from int to polynomials
    f2_int := orbit["f2int"];
    f3_int := orbit["f3int"];
    f2, f3 := convertToPolys(f2_int, f3_int, R, G, Bit2, Bit3);  // Now f2, f3 are in R

    // Permute x1, x2, x3 to get it in the correct form 
    if f2_int eq 1072 then 
        perm := [ x[2], x[1], x[3], x[4], x[5]];
        f2 := Evaluate(f2, perm);
        f3 := Evaluate(f3, perm);
    elif f2_int eq 656 then 
        perm := [ x[2], x[3], x[1], x[4], x[5]];
        f2 := Evaluate(f2, perm);
        f3 := Evaluate(f3, perm); 
    end if;

    orbit["f2int"] := Seqint([Integers() !Bit2(f2)[i] : i in [1..#Basis(V2)]], 2);
    orbit["f3int"] := Seqint([Integers() !Bit3(f3)[i] : i in [1..#Basis(V3)]], 2);
    new_orbits[key] := orbit;

    count +:= 1;
    if count mod 10000 eq 0 then count; end if;
    
end for;

containers := new_orbits;

// Save previous calculations
F := Open("orbits_lines.m", "w");
WriteObject(F, orbs_with_lines);
delete F;

F := Open("orbits_lines_contain.m", "w");
WriteObject(F, containers);
delete F;

F := Open("orbits_lines_secrtpt.m", "w");
WriteObject(F, rat_pters);
delete F;

F := Open("orbits_lines_sec.m", "w");
WriteObject(F, leftovers);
delete F;