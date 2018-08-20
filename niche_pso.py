import matplotlib.pyplot as plt 
import numpy as np
from benchmark_functions import *

step = 0.01
min_theta = -5.2
max_theta = 5.2
min_x = -3.4
max_x = 3.4



theta = np.arange(min_theta, max_theta, step)
x_s = np.arange(min_x, max_x, 0.01) 


for th_ in theta:
	values = [np.polyval([1, 0, -1, 0, 1, 0, th_, 0], x_) for x_ in x_s]
	#values = [np.polyval([-1, 0, 1, 0, th_, 0], x_) for x_ in x_s]
	#values = [np.polyval([-1, 0, 4, th_], x_) for x_ in x_s]
	#values = [np.polyval([-1, 1, th_], x_) for x_ in x_s]

	#values = [function4(x_, th_) for x_ in x_s]
	#values = [function2(x_, th_) for x_ in x_s]
	points = np.array([values, x_s])
	points_t = np.transpose(points)
	plt.grid()
	plt.plot(values, x_s, "-b")
	plt.xlim(min_x, max_x); plt.ylim(min_x, max_x)
	plt.xlabel(r'$\theta axis$'); plt.ylabel("X axis");
	plt.title("f(x, theta) = x^7 -x^5 + x^3 + x theta")
	#plt.title("f(x, theta) = -x^5 + x^3 + x theta")
	#plt.title(r'$f(x, \theta) = 4x-x^3+\theta$')
	#plt.title(r'$f(x, \theta) = x + \theta - x^2$')
	plt.show()
