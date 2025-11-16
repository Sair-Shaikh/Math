

ReadPointCountsCSV := function(file_paths)
    
    point_counts := AssociativeArray();

    for file_path in file_paths do 
        // Open the file
        F := Open(file_path, "r");

        // Loop through every line in the file
        while true do
            line := Gets(F);
            if IsEof(line) then break; end if;

            // Skip empty lines
            if line eq "" then
                continue;
            end if;

            // Split on commas
            fields := Split(line, ",");

            // First entry is the key
            key := Trim(fields[1]);

            // Convert the remaining entries into integers
            nums := [];
            for i in [2..#fields] do
                Append(~nums, StringToInteger(Trim(fields[i])));
            end for;

            // Store into the associative array
            point_counts[key] := nums;
            // printf "%o : %o\n", key, nums;
        end while;

        // Close the file and return
        delete F;
    end for;
    return point_counts;
end function;

file_paths := ["Dataset/point_counts_cubic.csv", "Dataset/point_counts_quad.csv"];
point_counts := ReadPointCountsCSV(file_paths);
orbit_paths := ["Dataset/orbits_lines_contain.m", "Dataset/orbits_lines_secrtpt2.m", "Dataset/orbits_lines_sec.m"];

printf "Number of orbits with point counts: %o\n", #point_counts;

new_orbits := AssociativeArray();
for file_path in orbit_paths do 
    F := Open(file_path, "r");
    orbits := ReadObject(F);
    delete F;

    for key in Keys(orbits) do 
        orbit := orbits[key];

        if not IsDefined(point_counts, key) then 
            key;
        else 
            orbit["pt_count"] := point_counts[key];
            new_orbits[key] := orbit;
        end if;
    end for; 
end for;


F := Open("Dataset/orbits_point_counts.m", "w");
WriteObject(F, new_orbits);
delete F;



