from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def recursive_states (mtrx, n_states):
	matrix_conect_states = np.zeros((n_states, n_states))

	# get connections between states
	P = mtrx
	for _ in range(n_states * 3):
		for i in range(n_states):
			matrix_conect_states[i] = np.copy([1. if P[i, j] > 0. else 0. for j in range(n_states)])
		P = np.dot(P, mtrx)

	return matrix_conect_states


def plot_niche_pso(positions_particles, main_swarm, th_, theta_s, x_s):
	positions_particles_t = np.transpose(positions_particles)

	for sb_w in main_swarm:
		# plot original graph
		plt.plot(theta_s, x_s, "-b")
		plt.grid()
		sb_particles = []
		if len(sb_w) > 0:
			for i in sb_w:
				sb_particles.append(positions_particles[i])
			sb_particles_t = np.transpose(sb_particles)
			plt.title(("Theta = %s" % (th_)), fontsize=16)

			plt.scatter(positions_particles_t[0], positions_particles_t[1], color='black', s=3.)
			plt.scatter(sb_particles_t[0], sb_particles_t[1], color='red', s=4.)
			#plt.xlim(-0.4, 0.4); plt.ylim(-2, 2);
			plt.show()


def plot_niche_pso2(positions_particles, th_, theta_s, x_s, reserved_nitches):

	plt.title(("Theta = %s" % (th_)), fontsize=16)
	plt.grid()
	if reserved_nitches != []:
		reserved_nitches_t = np.transpose(reserved_nitches)
		plt.scatter(reserved_nitches_t[0], reserved_nitches_t[1], color='green', s=8.)

	if positions_particles is not False:
		positions_particles_t = np.transpose(positions_particles)
		plt.scatter(positions_particles_t[0], positions_particles_t[1], color='red', s=2.)
	plt.plot(theta_s, x_s, "-b")
	#plt.xlim(-0.4, 0.4); plt.ylim(-2, 2);
	plt.show()


def plt_niche_pso3(mtrx_reserved_niches):
	lst_theta = []
	lst_ptos = []
	for i in range(len(mtrx_reserved_niches)):
		lst_theta.append(mtrx_reserved_niches[i][1])
		lst_ptos.append(mtrx_reserved_niches[i][2])

	lst_X = []
	lst_Y = []
	for i in range(len(lst_ptos)):
		aux = list(lst_ptos[i])

		for j in range(len(aux)):
			lst_Y.append(aux[j][1])
			lst_X.append(lst_theta[i])

	plt.grid()
	plt.scatter(lst_X, lst_Y, s=0.8, color="black")
	#plt.xlim(-0.4, 0.4); plt.ylim(-1.5, 1.5);
	plt.title("Niche PSO Algorithm", fontsize=16.)
	plt.xlabel('Theta axis'); plt.ylabel('X axis');
	plt.show()


