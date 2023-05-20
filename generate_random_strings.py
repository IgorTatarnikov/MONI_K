import random
import numpy as np


alphabet = ["A", "C", "G", "T"]
max_length = 20 
min_length = 8
strings = ["TAGATTACATTA"]
num_strings = 200

for i in range(num_strings):
    new_string = ""
    string_length = random.randint(min_length, max_length+1)

    for _ in range(string_length):
        new_string += random.choice(alphabet)

    strings.append(new_string)

f = open("test_strings.txt", "w")

for string in strings:
    f.write(f"{string}\n")

f.close()
