#this is a test just to make sure that the lc algorithm matches the lc alg as calculated by chris 

import l1_lc_alg as l1
import numpy as np
import math 
import matplotlib.pyplot as plt

d = np.linspace(0, 400, 1000)
p_max = []
p_min = []
for i in d:
    pos_max, pos_min, neg_max, neg_min = l1.light_cone(i)
    p_max.append(pos_max)
    p_min.append(pos_min)

plt.plot(d, p_max)
plt.plot(d, p_min)
plt.savefig("test.png")
