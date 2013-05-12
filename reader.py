lines = [line.split() for line in file("b.ab").readlines() if line.split()[0] == "1"]
print lines[0]
PAR = [float(line[4]) for line in lines]
L2 = [float(line[5]) for line in lines]
U1 = [float(line[6]) for line in lines]
U2 = [float(line[7]) for line in lines]
print U2
import matplotlib.pyplot as plt

plt.plot(PAR)
plt.show()

