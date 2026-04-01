
import random
import string

CHAPTERS = [1, 2, 3, 4]

EXERCISES_PER_SECTION = {
    1: [],
    2: [1, 10, 10, 16, 5], #, 10, 7],
    3: [2, 20, 0, 11, 6, 24, 7],
    4: [7, 0, 7, 5, 17],
}

DONE = {
    1: {},
    2: {"2.1.A", "2.1.B", 
        "2.2.E", "2.2.F", "2.2.J", 
        "2.3.A", "2.3.B", "2.3.C", "2.3.E", "2.3.G", "2.3.H", 
        "2.4.A", "2.4.B", "2.4.C", "2.4.D", "2.4.G", "2.4.H", "2.4.I", "2.4.J", "2.4.K"},
    3: {"3.2.D", "3.2.G", "3.2.M",
        "3.4.D", "3.4.J",
        "3.5.A", "3.5.B", "3.5.C", "3.5.D", "3.5.E"
        },
    4: {"4.1.A"},
}

UNIMPORTANT = {
    "2.3.D", 
    "3.2.B", "3.2.E", "3.2.H", "3.2.I", "3.2.R", "3.2.T",
    "3.4.A", "3.4.B", "3.4.C", "3.4.G"
    "4.1.B", "4.1.C"
}


ALL = set()

for c in CHAPTERS:
    for s, n in enumerate(EXERCISES_PER_SECTION.get(c, []), start=1):
        for letter in string.ascii_uppercase[:n]:
            ALL.add(f"{c}.{s}.{letter}")
ALL -= UNIMPORTANT

ALL_DONE = set([val for vals in DONE.values() for val in vals])

REMAINING = ALL - ALL_DONE - UNIMPORTANT

print(f"Total Questions: {len(ALL)}")
print(f"Questions Done: {len(ALL_DONE)}")
print(f"Questions left: {len(REMAINING)}")


print(f"\nRandom Choice: {random.choice(list(REMAINING))}\n")