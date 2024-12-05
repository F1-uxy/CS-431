# Better diff tool with tolerance

import sys

file1 = sys.argv[1]
file2 = sys.argv[2]
count = 0
with open(file1, 'r') as f1, open(file2, 'r') as f2:
    line_num = 1
    while True:
        l1 = f1.readline()
        l2 = f2.readline()

        if not l1 or not l2:
            break

        v1 = float(l1.strip())
        v2 = float(l2.strip())

        diff = abs(v1 - v2)

        if(diff > 1e-6):
            print(f"Line {line_num}: {v1} - {v2}")
            count += 1

    line_num += 1

print(f"Number of differences: {count}")
            
        
