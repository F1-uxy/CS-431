import matplotlib.pyplot as plt

times = [[],[],[]]
calculations = [[],[],[]]

with open('/home/cflux/Documents/Uni/Yr3/Intro to Parallel Computing/Homework/uoregon-cs431531-f24/homework01/results.txt', 'r') as file:
    for line in file:
        pairs = line.strip().split()
        for idx, pair in enumerate(pairs):
            time, calc = map(float, pair.split(':'))
            times[idx % 3].append(time)
            calculations[idx % 3].append(calc)


plt.boxplot(times[1])
plt.title("Execution Time for Pi Calculation")
plt.ylabel("Time (seconds)")
plt.show()
