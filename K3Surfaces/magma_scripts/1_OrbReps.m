

/*
    Use Magma's built=in for computing orbit representatives and sizes
*/
function GetOrbitRepresentatives(V)
    GV := ActionGroup(V);
    phi := GModuleAction(V);
    reps := [<o[1], #Orbit(GV, o[1])> : o in Orbits(GV)];
    return reps;   
end function;


/*
    Return 1 if K3 surface is smooth, 0 otherwise.
*/
function IsSmooth(P, R, f, g)
    S := Scheme(P, [f, g]);
    if IsNonsingular(S) then
        return 1;
    else 
        return 0;
    end if;
end function;


/*
    Use Burnside's Lemma to enumerate expected number of orbits -- to check against
*/

function BurnsideOrbCount(V)
    count := 0;
    conjugacy_classes := ConjugacyClasses(ActionGroup(V));
    for i in [1..#conjugacy_classes] do
        cls := conjugacy_classes[i];
        length := cls[2];
        rep := cls[3];
        Mat_rep := Matrix(rep);
        nullity := Dimension(Kernel(Mat_rep-IdentityMatrix(BaseRing(V), Dimension(V))));
        X_rep := 2^nullity;
        count := count+length*X_rep;
    end for;
    
    return count/#ActionGroup(V);
end function;


function OrbitsByFiltration(V, F)
    
    // Setup Array:
    OrbReps := [];
    
    if IsEmpty(F) then
        return GetOrbitRepresentatives(V);
    end if;

    "Level ", #F;
    
    // Take the last factor in the filtration 
    B := F[#F];  
    
    // Calculate the orbits of quotient V/B
    U, q_VmodB := quo<V |B>;
    OrbRepsU := OrbitsByFiltration(U, Prune(F));

    iotaBV := Morphism(B, V);
    phi := GModuleAction(U);
    Gsize := #ActionGroup(V);
    
    // For each quotient orbit, get orbits of V
    for j in [1..#OrbRepsU] do
        "    Orbit ", j, " of ", #OrbRepsU;
        y := OrbRepsU[j][1];
        x :=  y @@ q_VmodB;
        fiber_over_y := {x + iotaBV(b) : b in B}; // Preimage is disj. union. of orbits
        
        yStab := Stabilizer(ActionGroup(U), y) @@ phi; // The stabilizer as a subgroup of G
        
        if #yStab eq 1 then 
            for p in fiber_over_y do
                Append(~OrbReps, <p, Gsize>);
            end for;
            continue;
        end if;
        
        V_yStab := Restriction(V, yStab); // Restrict action on V to stabv        
        GVstab := ActionGroup(V_yStab);

        A := AssociativeArray();
        for u in fiber_over_y do
            A[u] := false;
        end for;
    
        for u in fiber_over_y do
            if A[u] eq true then continue; end if;
            orbit := Orbit(GVstab, u);
            rep := (V!u);
            
            Append(~OrbReps, <rep,  Gsize/(#GVstab/#orbit)>);
            
            // Update A to mark progress
            for a in orbit do
                A[a] := true;
            end for;
        end for;
    end for; 

    return OrbReps;

end function;

// Setup
k := FiniteField(2);
G := GL(5, k);
R<[x]> := PolynomialRing(k, 5);
P := ProjectiveSpace(R);
V2, Bit2 := GModule(G, R, 2);

// Get all non-zero orbits of V2
V2_orbreps := GetOrbitRepresentatives(V2);
V2_orbreps := [o : o in V2_orbreps | o[1] ne 0 ];

// Initialize array to store all sextic K3 orbreps
K3_orbreps := AssociativeArray();
K3_smooth := AssociativeArray();

smooth_count := 0;
bs_count := 0;
for i in [1..#V2_orbreps] do
    // Setup
    curr_count := 0;
    "V2 Orbit #", i, "out of", #V2_orbreps;
    current_orbit := AssociativeArray();
    
    f := V2!V2_orbreps[i][1];
    
    // Create V3/fL acted on by fStab
    phi := GModuleAction(V2);
    fStab := (Stabilizer(ActionGroup(V2), f)) @@ phi;
    
    V3, Bit3 := GModule(fStab, R, 3);   // All deg3 hom. poly.
    L, BitL := GModule(fStab, R, 1);    // All deg1 hom. poly.
    fL := sub<V3| [Bit3((f@@Bit2)*(l @@ BitL)) : l in Basis(L)]>; 
    V3modfL, q_V3modfL := quo<V3 | fL>;

    // Count orbits using Burnside
    curr_total := BurnsideOrbCount(V3modfL);
    bs_count +:= curr_total;

    // Create one step filtration for V3
    subs := Submodules(V3modfL);
    B := subs[Floor(#subs/2)];
    "Submodule: ", B;

    // Find orbits in V3modFL
    V3modfL_orbits := OrbitsByFiltration(V3modfL, [B]);
    "Expected Orbits: ", curr_total;
    "Total Orbits Found: ", #V3modfL_orbits;
    for j in [1..#V3modfL_orbits] do
        orb := V3modfL_orbits[j];
        rep := orb[1]; 
        size := orb[2];
        g := rep @@ q_V3modfL;

        idx := Sprint(i) cat ("-" cat Sprint(j));
        fint := Seqint([Integers() !f[i] : i in [1..#Basis(V2)]], 2);
        gint := Seqint([Integers() !g[i] : i in [1..#Basis(V3)]], 2);
        smooth := IsSmooth(P, R, f @@ Bit2, g @@ Bit3);
        current_orbit["idx"] := idx;
        current_orbit["f2int"] := fint;
        current_orbit["f3int"] := gint;
        current_orbit["smooth"] := smooth;
        current_orbit["size"] := size;
        K3_orbreps[idx] := current_orbit;

        if smooth eq 1 then 
            K3_smooth[idx] := current_orbit;
            smooth_count +:= 1;
        end if;

        curr_count +:= 1;
        if curr_count mod 50000 eq 0 then
            "Done:", curr_count;
        end if;

        if smooth_count mod 10000 eq 0 and smooth_count ge 1 then
            "Smooth:", smooth_count;
        end if;
    
    end for;


    // Save previous calculations
    F := Open("orbit_reps_all.m", "w");
    WriteObject(F, K3_orbreps);
    delete F;

    F := Open("orbit_reps_smooth.m", "w");
    WriteObject(F, K3_smooth);
    delete F;
end for;




