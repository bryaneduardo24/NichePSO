import numpy as np
import random

#dx/dt = 4x-x^3+c where 0 > = c <= 4
"""
min_theta = -4
max_theta = 4
min_x = -4
max_x = 4
"""
def function1(x, c):
	# [-1, 0, 4, c]
	return (4 * x) - (x**3) + c

#dx/dt = x(1 - x) - c where  -1 > c < 1
"""
min_theta = -2
max_theta = 2
min_x = -4
max_x = 4

"""
def function2(x, c):
	# [-1, 1, -c]
	return x * (1 - x) - c

#dx/dt = x^2 - cx + 4
def function3(x, c):
	# [1, -c, 4]
	return x**2 - (c * x) + 4
	
#dx/dt = x + c - x**2
"""
step = 0.1
min_theta = -10
max_theta = 10
min_x = -10
max_x = 10
"""
def function4(x, c):
	# [-1, 1, c]
	return x + c - x**2

#f(x) = 2x^4 - x^3 - 2x^2 
#def function5(x):
#	return 2 * (x**4) - x**3 - 2 * (x**2)

#f(x) = x^6 - 3.2x^4 + 2x^2
#def function6(x):
#	return x**6 - 3.2 * (x**4) + 2 (x**2)

#f(x) = x^4 - 2x^2 
#def function7(x):
#	return x**4 - 2 * (x**2)
	
#cross in tray modified function where -10 > x < 10	
# FIXME : Que puedo variar?
def crossintray(x):
	return -0.0001 * (np.sin(x) * np.sin(x)) * (np.exp(abs(100 - np.sqrt(x**2 + x**2)/np.pi)) + 1)**2

#eggholder modified function where -80 > x < 20
# FIXME : Que puedo variar?
def eggholder(x):
	return - (x + 47) * np.sin(np.sqrt(abs(x + x/2 + 47)) - x * np.sin(np.sqrt(abs(x - (x + 47)))))
	
