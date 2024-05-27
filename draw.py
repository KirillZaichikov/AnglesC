import numpy as np
import matplotlib.pyplot as plt

size=5
eps=0.0001
eps1=0.00001
f = open("build/File.txt", "r")
fig, ax = plt.subplots()
ax.plot()
xChaosP = []
yChaosP = []
xChaosQ = []
yChaosQ = []
xStable = []
yStable = []
xCircle = []
yCircle = []
mas = f.read()
mas = mas.split(", ")
mas = mas[:-1]
for a in mas:
    a = a.split(" ")
    print(a)
    if float(a[2])>eps:
        if float(a[3])>eps1:
            xChaosP.append(float(a[0]))
            yChaosP.append(float(a[1]))
        else:
            xChaosQ.append(float(a[0]))
            yChaosQ.append(float(a[1]))
    elif float(a[2])<-eps:
        xStable.append(float(a[0]))
        yStable.append(float(a[1]))
    else:
        xCircle.append(float(a[0]))
        yCircle.append(float(a[1]))
print(xStable, yStable)
ax.scatter(xChaosP, yChaosP, c='orange', s=size)
ax.scatter(xChaosQ, yChaosQ, c='blue', s=size)
ax.scatter(xStable, yStable, c='#42aaff', s=size)
ax.scatter(xCircle, yCircle, c='white', s=size)
plt.show()

