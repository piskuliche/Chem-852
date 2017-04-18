import numpy as np

for t in range(0,100000):
    time=t/100.0
    x=np.cos(.2*time)
    print("%s %s" % (time, x))
