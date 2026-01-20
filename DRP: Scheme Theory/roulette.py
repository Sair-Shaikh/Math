
import random

CHAPTERS = [1, 2, 3, 4]

EXERCISES_PER_SECTION = {
    1: [],
    2: [1, 10, 10, 16, 5, 10, 7],
    3: [],
    4: [],
}

DONE = {
    1: [],
    2: {"2.1.A", "2.2.F", "2.2.J"},
    3: [],
    4: [],
}



while True:
    c = random.choice(CHAPTERS)
    if not len(EXERCISES_PER_SECTION[c]): continue

    s = random.randint(1, len(EXERCISES_PER_SECTION[c]))

    e = random.randint(1, EXERCISES_PER_SECTION[c][s - 1])

    letter = chr(64+e)

    exercise = f"{c}.{s}.{letter}"

    if exercise in DONE[c]: 
        continue
    else: 
        print(f"{c}.{s}.{letter}")
        break