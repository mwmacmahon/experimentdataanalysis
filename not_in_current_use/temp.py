# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 11:29:31 2016

@author: mwmac
"""

# %%

import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-2, 2, 100)
x1 = x[np.argwhere(x < -0.2)]
x2 = x[np.argwhere(np.logical_and(x >= -0.2, x < 0.2))]
x3 = x[np.argwhere(x >= 0.2)]

y1 = np.clip(-1.5/(abs(x1) + 0.12), -3, 0)
y2 = np.interp(x2, (-0.2, 0.2), (-3, 1))
y3 = np.clip(1.1/(abs(x3) + 0.8), 0, 1)
y = np.concatenate([y1, y2, y3])

plt.plot(x, y, 'd')