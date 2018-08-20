#------------------------------------------------------------------------------+
# DEP-FIE
# Author : Bryan E. Martinez
# July, 2017
#
# This script was created using a code for particles of:
#   Nathan A. Rooy
#   Simple Particle Swarm Optimization (PSO) with Python
#   July, 2016
#
#------------------------------------------------------------------------------+
from __future__ import division
import random
import math
import numpy as np
from scipy.spatial import distance
from tools import *
import time


###############################################################################
# F => Mutation in range (0,2), CR => Recombination in range (0,1)
def differential_evolution(bounds, x, v, fitness_function, th, d, f=.4,cr=.1,maxiter=10):
    it=0
    err_best_g=-1
    pos_best_g=[]
    len_x=len(x)
    while stop_condition(it,maxiter):
        v=mutation(len_x,x,v,f,d,fitness_function,th)
        v=recombination(len_x,x,v,cr,d,fitness_function,th)

        for j in range(0, len_x):
            x[j].evaluate(fitness_function, th)

            # determine if current particle is the best (globably)
            if x[j].err_i<err_best_g or err_best_g==-1:
                err_best_g=float(x[j].err_i)
                pos_best_g=list(x[j].position_i)

        x=selection(len_x,x,v)

        it+=1


    for j in range(0, num_particles):
        lst_particles[j].update_position(bounds)
    return (err_best_g, pos_best_g)

def stop_condition(it,maxiter):
    if maxiter!=0:
        if it<maxiter:
            return True
        else:
            return False

def mutation(len_x,x,v,f,d,fitness,th):
    for i in range(len_x):
        r=[]
        for j in range(3):
            add_r(len_x,r)
        for j in range(d):
            v[i].position_i[j]=float(x[r[0]].position_i[j])+f*(float(x[r[1]].position_i[j]) - float(x[r[2]].position_i[j]))
            v[j].evaluate(fitness, th)
    return v


def recombination(len_x,x,v,cr,d,fitness,th):
#    import ipdb; ipdb.set_trace()
    for i in range(len_x):
#            j=int(math.floor(random.uniform(0,d)))
        j=1
        if random.uniform(0,1)<=cr:
            v[i].position_i[j]=x[i].position_i[j]
            v[j].evaluate(fitness, th)
    
    return v


def selection(len_x,x,u):
    for i in range(len_x):
        if u[i].err_i<x[i].err_i:
            x[i]=u[i]
    
    return x

#Obtain indices r1, r2 and r3 such that are distinct
def add_r(len_x,r):
    tmp_r=int(math.floor(random.uniform(0,len_x)))
    if tmp_r not in r:
        r.append(tmp_r)
    else:
        add_r(len_x,r)


###############################################################################


def get_ms(lst_points, points=2):
    n = len(lst_points)

    lst_ = [np.array(lst_points[i + 1]) - np.array(lst_points[i]) for i in range(n - 1)]

    if points == 2:
        lst_ = list(lst_[-1])
        if lst_[0] == 0:
            m = lst_[1]
        else:
            m = lst_[1] / lst_[0]
    else:
        m = -1000000000000
    return m


def print_mtrx(mtrx_reserved_niches):
    print "\t"
    print "[[num roots0, value theta0, niches0, error niches0], [], ...]"
    for i in range(len(mtrx_reserved_niches)):
        print ("%s - %s" % (i, mtrx_reserved_niches[i]))
    print "<=====================>"


def explore_with_bisection(mtrx_reserved_niches):
    is_change = False
    pos_critical = []
    for i in range(len(mtrx_reserved_niches)):
        if i > 0:
            auxp0 = mtrx_reserved_niches[i - 1]
            auxp1 = mtrx_reserved_niches[i]

            if auxp0[0] != auxp1[0]:
                pos_critical.append([i - 1, i])

    for k in range(len(pos_critical) - 1, -1, -1):
        i, j = pos_critical[k]

        p0 = mtrx_reserved_niches[i]
        p1 = mtrx_reserved_niches[j]

        lefts = np.copy(p0[0])
        rigths = np.copy(p1[0])
        num_roots = list([p0[0], p1[0]])
        maxroots = np.max(num_roots)
        minroots = np.min(num_roots)

        diff_between_ps = p1[1] - p0[1]
        while diff_between_ps >= th_bisection:
            # get middle point to evaluate new theta value
            middle = (diff_between_ps / 2.) + p0[1]

            # get theta values
            theta_s = [function_to_val(x_, middle) for x_ in x_s]

            # first run NichePSO Algorithm ...
            lst_pos_particles = create_particles_position([(0, 0), (min_x, max_x)], type_bound=2, n_rows=n_particles_niche, n_cols=0.)

            reserved_niches_err, reserved_niches = NichePSO(fitness, lst_pos_particles, [(0, 0), (min_x, max_x)], maxiter, middle, theta_s, x_s)

            new_roots_f = len(reserved_niches)

            # print ("There are %s roots range, and now get %s" % ([minroots, maxroots], new_roots_f))
            if new_roots_f < minroots or new_roots_f > maxroots:
                diff_between_ps = 0.0

            else:
                if new_roots_f == rigths:
                    # update range right to left
                    if j != len(pos_critical) - 1:
                        mtrx_reserved_niches[j] = list([new_roots_f, middle, reserved_niches, reserved_niches_err])
                    else:
                        mtrx_reserved_niches.append(list(mtrx_reserved_niches[-1]))
                        mtrx_reserved_niches[-2] = list([new_roots_f, middle, reserved_niches, reserved_niches_err])
                        is_change = True
                        break

                elif new_roots_f == lefts:
                    # update range left to rigth
                    if i != 0:
                        mtrx_reserved_niches[i] = list([new_roots_f, middle, reserved_niches, reserved_niches_err])
                    else:
                        mtrx_reserved_niches.insert(0, list(mtrx_reserved_niches[0]))
                        mtrx_reserved_niches[1] = list([new_roots_f, middle, reserved_niches, reserved_niches_err])
                        is_change = True
                        break

                else:
                    # append new route
                    mtrx_reserved_niches.insert(i + 1, list(mtrx_reserved_niches[i]))
                    mtrx_reserved_niches.insert(i + 2, list([new_roots_f, middle, reserved_niches, reserved_niches_err]))
                    mtrx_reserved_niches.insert(i + 3, list(mtrx_reserved_niches[j + 2]))
                    is_change = True
                    break

                # update
                p0 = mtrx_reserved_niches[i]
                p1 = mtrx_reserved_niches[j]
                diff_between_ps = p1[1] - p0[1]
    return is_change, mtrx_reserved_niches


def check_explore_with_bisection(mtrx_reserved_niches):
    #print_mtrx(mtrx_reserved_niches)

    is_change = True
    while is_change:
        #print_mtrx(mtrx_reserved_niches)
        is_change, mtrx_reserved_niches = explore_with_bisection(mtrx_reserved_niches)
        #print_mtrx(mtrx_reserved_niches)
    return mtrx_reserved_niches


def create_particles_position(bounds, type_bound=1, n_rows=0., n_cols=0., threshold=1., m_point=0.):
    positions_particles = []
    x_range = []
    y_range = []

    # [(min_theta2, max_theta2), (min_x, max_x)]
    if type(bounds) == list and type_bound == 2:
        x1_bound = bounds[0]
        x2_bound = bounds[1]

        dst_x = x1_bound[1] - x1_bound[0]
        dst_y = x2_bound[1] - x2_bound[0]

        if n_cols != 0.:
            x_step = dst_x / float(n_cols)
            x_range = np.arange(x1_bound[0], x1_bound[1] + x_step, x_step)
        else:
            x_range = np.array([0])

        if n_rows != 0.:
            y_step = dst_y / float(n_rows)
            y_range = np.arange(x2_bound[0], x2_bound[1] + y_step, y_step)
        else:
            y_range = np.array([0])

    elif type(bounds) == list and type_bound == 1:
        x1_bound = bounds

        if n_cols != 0.:
            pi = x1_bound[0] - threshold
            pf = x1_bound[0] + threshold
            x_range = np.array([random.uniform(pi, pf) for _ in range(int(n_cols))])
        else:
            x_range = np.array([0])

        if n_rows != 0.:
            pi = x1_bound[1] - threshold
            pf = x1_bound[1] + threshold
            y_range = np.array([random.uniform(pi, pf) for _ in range(int(n_rows))])
        else:
            y_range = np.array([0])

    elif type(bounds) == list and type_bound == 0:
        x1_bound = bounds

        if n_cols != 0.:
            # FIXME : needed limits for a distribution 
            pi = x1_bound[0] - threshold
            pf = x1_bound[0] + threshold
            x_range = np.array([random.uniform(pi, pf) for _ in range(int(n_cols))])
        else:
            x_range = np.array([0])

        if n_rows != 0.:
            p_ = x1_bound[1]
            y_range = [p_]

            # for 3 positions distributed we need 5 rows ...

            if m_point >= 0:
                y_range.append(p_ - threshold)
                for i in range(1, int(n_rows) - 1):
                    y_range.append(p_ + threshold / i)

            else:
                y_range.append(p_ + threshold)
                for i in range(1, int(n_rows) - 1):
                    y_range.append(p_ - threshold / i)
            #y_range = np.array([random.uniform(pi, pf) for _ in range(int(n_rows))])
        else:
            y_range = np.array([0])


    for i in x_range:
        for j in y_range:
            positions_particles.append([i, j])

    # positions_particles_t = np.transpose(positions_particles)
    # plt.scatter(positions_particles_t[0], positions_particles_t[1], s=3.0)
    # plt.show()

    return positions_particles


def threshold_niche(point_niche, new_point, rad_th):
    # calculate euclidean distace
    x_hat = new_point[0] - point_niche[0]
    y_hat = new_point[1] - point_niche[1]
    dst_ = np.sqrt( x_hat ** 2 + y_hat ** 2)
    if dst_ < rad_th:
        theta = np.arctan2(y_hat, x_hat)
        x_new = np.cos(theta) * rad_th + point_niche[0]
        y_new = np.sin(theta) * rad_th + point_niche[1]

    else:
        x_new = new_point[0]
        y_new = new_point[1]

    return [x_new, y_new]


def get_pos_particles_from_lst(lst_particles):
    num_particles = len(lst_particles)
    positions_p = []
    for j in range(0,num_particles):
        positions_p.append(lst_particles[j].position_i)
    return positions_p


def get_particles(lst_particles, lst_positions):
    aux_positions = [lst_particles[pos] for pos in lst_positions]
    return aux_positions


def add_new_sub_swarms(lst_main_swarm, positions_particles, lst_err_best_g, lst_pos_best_g):
    if len(lst_main_swarm[-1]) >= 2:
        rename_particles = np.copy(lst_main_swarm[-1])
        lst_main_swarm = np.delete(lst_main_swarm, -1)
        lst_err_best_g.pop()
        lst_pos_best_g.pop()

        new_sub, new_lst_err_best_g, new_lst_pos_best_g = create_sub_swarms(positions_particles, rename_particles)

        lst_main_swarm = np.concatenate((lst_main_swarm, new_sub), axis=0)

        for i, j in zip(new_lst_err_best_g, new_lst_pos_best_g):
            lst_err_best_g.append(i)
            lst_pos_best_g.append(j)

    return (lst_main_swarm, lst_err_best_g, lst_pos_best_g)


def create_sub_swarms(positions_particles, rename_particles=False):
    n = len(positions_particles)
    mtrx_dst_particles = np.zeros((n,n))
    for i, particle in enumerate(positions_particles):
        for j in range(i + 1, n):
            particle2 = positions_particles[j]
            dst_between_particles = distance.euclidean(particle, particle2)
            # create a link in mtrx if dst between particles is le to dst_max_particles
            if dst_between_particles <= dst_max_particles:
                mtrx_dst_particles[i][j] = 1.  # dst_between_particles
                mtrx_dst_particles[j][i] = 1.  # dst_between_particles

    # get connections between states
    mtrx_dst_particles = recursive_states(mtrx_dst_particles, n)

    # get swarms
    main_swarm = []
    without_sub_swarm = []
    is_used_particle = [0 for i in range(n)]
    for i in range(n):
        aux_vtr = np.copy(mtrx_dst_particles[i])
        aux_vtr[i] = 0.
        aux_sum_vtr = np.sum(aux_vtr)

        if aux_sum_vtr == 0:
            if is_used_particle[i] == 0:
                without_sub_swarm.append(i)

        elif aux_sum_vtr > 0:
            sub_swarm = []
            sub_swarm.append(i)
            is_used_particle[i] = 1

            for j in range(i + 1, n):
                if mtrx_dst_particles[i][j] == 1.:
                    sub_swarm.append(j)
                    is_used_particle[j] = 1

            for k in range(n):
                for l in sub_swarm:
                    mtrx_dst_particles[k][l] = 0.

            main_swarm.append(sub_swarm)

    main_swarm.insert(0, [])
    main_swarm.append(without_sub_swarm)

    # Return a numpy matrix ...
    main_swarm = [np.array(s_w, dtype=int) for s_w in main_swarm]
    main_swarm = np.array(main_swarm)

    # mapping is true
    if type(rename_particles) is np.ndarray and len(rename_particles) > 0:
        for i in range(len(main_swarm)):
            set_p = np.copy(main_swarm[i])
            for j in range(len(set_p)):
                n_pos = set_p[j]
                main_swarm[i][j] = int(rename_particles[n_pos])

    # return other vectors
    rows = len(main_swarm)
    lst_err_best_g = [-1 for i in range(rows)]
    lst_pos_best_g = [[] for i in range(rows)]

    return (main_swarm, lst_err_best_g, lst_pos_best_g)


def send_lonely_particles_to_main_swarm(lst_main_swarm, lst_err_best_g, lst_pos_best_g):
    n = len(lst_main_swarm)
    lst_aux = []

    for i in range(1, n - 1):
        sub_swarm = lst_main_swarm[i]
        if len(sub_swarm) == 1:
            lst_aux.append(i)

    for i in range(len(lst_aux) - 1, -1, -1):
        idx = lst_aux[i]
        lst_main_swarm[-1] = np.concatenate((lst_main_swarm[-1], lst_main_swarm[idx]), axis=0)
        lst_main_swarm = np.delete(lst_main_swarm, idx)
        lst_err_best_g.pop(idx)
        lst_pos_best_g.pop(idx)

    return lst_main_swarm, lst_err_best_g, lst_pos_best_g


def clean_sub_swarms(lst_main_swarm, lst_err_best_g, lst_pos_best_g):
    n = len(lst_main_swarm)

    # remove the lst null
    lst_aux = []
    for i in range(1, n - 1):
        lst_ = lst_main_swarm[i]
        if len(lst_) == 0:
            lst_aux.append(i)

    for i in range(len(lst_aux) - 1, -1, -1):
        idx = lst_aux[i]
        lst_main_swarm = np.delete(lst_main_swarm, idx)
        lst_err_best_g.pop(idx)
        lst_pos_best_g.pop(idx)

    return lst_main_swarm, lst_err_best_g, lst_pos_best_g


# 1 -> without_sub_swarm to sub_swarm
# 2 -> sub_swarm to without_sub_swarm
def is_particle_sub_swarms(lst_main_swarm, lst_particles, opt):
    without_sub_swarm = lst_main_swarm[-1]
    lwss = len(without_sub_swarm)

    is_change = False
    for i in range(len(lst_main_swarm) - 1):
        sub_swarm = lst_main_swarm[i]
        if is_change == True:
            break

        for k in range(len(sub_swarm)):
            sub_particle = sub_swarm[k]

            for l in range(lwss):
                wout_particle = without_sub_swarm[l]
                p1 = lst_particles[sub_particle]
                p2 = lst_particles[wout_particle]
                dst_between_particles = distance.euclidean(p1, p2)

                if dst_between_particles <= dst_max_particles:
                    if opt == 1:
                        aux0_ = np.copy(sub_swarm)
                        aux0_ = np.append(aux0_, wout_particle)
                        lst_main_swarm[i] = np.copy(aux0_)

                        aux1_ = np.copy(without_sub_swarm)
                        aux1_ = np.delete(aux1_, l)
                        lst_main_swarm[-1] = np.copy(aux1_)

                        is_change = True
                        break

                    elif opt == 2:
                        aux1_ = np.copy(without_sub_swarm)
                        aux1_ = np.append(aux1_, sub_particle)
                        lst_main_swarm[-1] = np.copy(aux1_)

                        aux0_ = np.copy(sub_swarm)
                        aux0_ = np.delete(aux0_, k)
                        lst_main_swarm[i] = np.copy(aux0_)

                        is_change = True
                        break

    return (is_change, lst_main_swarm)


def check_particle_sub_swarms(lst_main_swarm, lst_particles, opt):
    is_change = True
    while is_change:
        is_change, lst_main_swarm = is_particle_sub_swarms(lst_main_swarm, lst_particles, opt)
    return lst_main_swarm


def merge_sub_swarms(lst_main_swarm, lst_particles):
    lst_main_swarm = np.array(lst_main_swarm)
    is_change = False
    n = len(lst_main_swarm) - 1

    for i in range(n):
        for j in  range(i + 1, n):
            sub_swarm = np.copy(lst_main_swarm[i])
            sub_swarm2 = np.copy(lst_main_swarm[j])

            for k in range(len(sub_swarm)):
                for l in range(len(sub_swarm2)):
                    p1 = lst_particles[sub_swarm[k]]
                    p2 = lst_particles[sub_swarm2[l]]
                    dst_between_particles = distance.euclidean(p1, p2)

                    if dst_between_particles <= dst_max_particles:
                        aux0_ = sub_swarm
                        aux1_ = sub_swarm2
                        aux0_ = np.concatenate((aux0_, aux1_), axis=0)

                        lst_main_swarm[i] = np.copy(aux0_)
                        lst_main_swarm[j] = np.copy([])

                        is_change = True
                        break

    return (is_change, lst_main_swarm)


def check_merge_sub_swarms(lst_main_swarm, lst_particles):
    is_change = True
    while is_change:
        is_change, lst_main_swarm = merge_sub_swarms(lst_main_swarm, lst_particles)
    return lst_main_swarm


# function we are attempting to optimize (minimize)
def fitness(point, th_):
    # [theta, x]
    coor0 = point[0]
    coor1 = point[1]

    total = np.abs(function_to_val(coor1, th_))
    # total = 1. / (1. + total)

    # print (("theta : %s, x : %s, th : %s, total : %s") % (coor0, coor1, th_, total))
    return total

"""
# function we are attempting to optimize (minimize)
def fitness(x, th_):
    total=0
    for i in range(len(x)):
        total+=x[i]**2
    return total
"""


def get_fuzzy_niches(lst_particles, fitness_function, th_):

    idx_best_error = []
    idx_best_positions = []

    llp = len(lst_particles)
    for i in range(llp):
        lst_particles[i].evaluate(fitness_function, th_)
        idx_best_error.append(lst_particles[i].err_i)
        idx_best_positions.append(lst_particles[i].position_i)

    pos_min = []
    for i in range(1, llp - 1):
        m0 = (idx_best_error[i] - idx_best_error[i - 1])
        m1 = (idx_best_error[i + 1] - idx_best_error[i])
        if m0 < 0 and m1 >= 0:
            pos_min.append(i)

    # create sub swarms
    lstms = [i for i in range(llp)]
    for i in range(len(pos_min) - 1, -1, -1):
        lstms.pop(pos_min[i])

    # create main
    lst_main_s = []
    rows = len(pos_min) + 2
    lst_err_best_g = [-1 for i in range(rows)]
    lst_pos_best_g = [[] for i in range(rows)]

    i = 1
    lst_main_s.append(np.array([], dtype=int))
    for val in pos_min:
        lst_main_s.append(np.array([val], dtype=int))
        lst_err_best_g[i] = idx_best_error[val].copy()
        lst_pos_best_g[i] = list(lst_particles[val].position_i)
        i += 1
    lst_main_s.append(np.array(lstms, dtype=int))
    lst_main_s = np.array(lst_main_s)

    return (lst_main_s, lst_err_best_g, lst_pos_best_g)


# FIXME : can be deleted ...
"""
def get_best_of_without_sub_swarm(lst_main_swarm, lst_particles, fitness_function, th_):
    lst_particles_ = get_particles(lst_particles, lst_main_swarm[-1])
    llp = len(lst_particles_)
    for i in range(llp):
        lst_particles_[i].evaluate(fitness_function, th_)

    idx_best_positions = []
    idx_best_error = []
    for i in range(llp - 1, -1, -1):
        if lst_particles_[i].err_i < 0.02:
            lst_main_swarm[-1] = np.delete(lst_main_swarm[-1], i)
            idx_best_positions.append(lst_particles_[i].position_i)
            idx_best_error.append(lst_particles_[i].err_i)

    return lst_main_swarm, idx_best_error, idx_best_positions
"""

def get_particles_nearest(lst_main_swarm, lst_pos, threshold_):
    for i in range(len(lst_main_swarm) - 1):
        sub_swarm = lst_main_swarm[i]
        if len(sub_swarm) >= 1:

            without_sub_swarm = lst_main_swarm[-1]
            lst_idx = []
            for j in range(len(without_sub_swarm)):
                p0 = lst_pos[sub_swarm[0]]
                p1 = lst_pos[without_sub_swarm[j]]
                new_dst = distance.euclidean(p0, p1)
                if new_dst <= threshold_:
                    lst_pos[without_sub_swarm[j]] = list(p0)
                    lst_idx.append(j)

            for j in range(len(lst_idx) - 1, -1, -1):
                idx = lst_idx[j]
                lst_main_swarm[i] = np.append(lst_main_swarm[i], without_sub_swarm[idx])
                lst_main_swarm[-1] = np.delete(lst_main_swarm[-1], idx)


    return lst_pos, lst_main_swarm


def remove_repeat(lst_err, lst_pos, step_theta_):
    n = len(lst_err)
    for i in range(n - 1):
        for j in range(i + 1, n):
            p1 = lst_pos[i]
            p2 = lst_pos[j]
            if len(p1) > 0 and len(p2) > 0:
                dst = distance.euclidean(p1, p2)
                if dst < step_theta_:
                    if lst_err[i] < lst_err[j]:
                        lst_err[i] = -1
                        lst_pos[i] = list([])
                    else:
                        lst_err[j] = -1
                        lst_pos[j] = list([])

    for i in range(n - 1, -1, -1):
        if lst_err[i] == -1:
            lst_err.pop(i)
            lst_pos.pop(i)

    return lst_err, lst_pos


def filter_to_minimal(lst_err_best_g, lst_pos_best_g, err_):
    reserved_niches_err = []
    reserved_niches = []
    for i in range(len(lst_err_best_g)):
        lebg = lst_err_best_g[i]
        if lebg <= err_:
            reserved_niches_err.append(lst_err_best_g[i])
            reserved_niches.append(lst_pos_best_g[i])

    return reserved_niches_err, reserved_niches


def PSO(num_particles, bounds, lst_point_niches, lst_particles, fitness_function, err_best_g, pos_best_g, th_, w=0.5, c1=1, c2=2):
    #
    # reject a particle near to a niche previous
    # FIXME : Not compare all vs all particles, only the nearest particles around one
    for i in range(0, len(lst_point_niches)):
        point_niche = lst_point_niches[i]
        for j in range(0, num_particles):
            new_point = lst_particles[j].position_i
            new_point = threshold_niche(point_niche, new_point, rad_th=reject_particles_error)
            lst_particles[j].position_i = new_point


    # cycle through particles in swarm and evaluate fitness
    for j in range(0, num_particles):
        lst_particles[j].evaluate(fitness_function, th_)

        # determine if current particle is the best (globably)
        if lst_particles[j].err_i<err_best_g or err_best_g==-1:
            err_best_g=float(lst_particles[j].err_i)
            pos_best_g=list(lst_particles[j].position_i)

    # cycle through swarm and update velocities and position
    for j in range(0, num_particles):
        lst_particles[j].update_velocity(pos_best_g, w, c1, c2)
        lst_particles[j].update_position(bounds)

    return (lst_particles, err_best_g, pos_best_g)


class Particle:
    def __init__(self,x0):
        self.position_i=[]          # particle position
        self.velocity_i=[]          # particle velocity
        self.pos_best_i=[]          # best position individual
        self.err_best_i=-1          # best error individual
        self.err_i=-1               # error individual

        for i in range(0,num_dimensions):
            self.velocity_i.append(random.uniform(-0.5, 0.5))
            self.position_i.append(x0[i])

    # evaluate current fitness
    def evaluate(self,costFunc, th_):
        self.err_i=costFunc(self.position_i, th_)

        # check to see if the current position is an individual best
        if self.err_i<self.err_best_i or self.err_best_i==-1:
            self.pos_best_i=self.position_i
            self.err_best_i=self.err_i

    # update new particle velocity
    def update_velocity(self,pos_best_g, w, c1, c2):
        # w  = 0.5      # constant inertia weight (how much to weigh the previous velocity)
        # c1 = 1        # cognative constant
        # c2 = 2        # social constant
        for i in range(0,num_dimensions):
            r1=random.random()
            r2=random.random()

            vel_cognitive=c1*r1*(self.pos_best_i[i]-self.position_i[i])
            vel_social=c2*r2*(pos_best_g[i]-self.position_i[i])
            self.velocity_i[i]=w*self.velocity_i[i]+vel_cognitive+vel_social

    # update the particle position based off new velocity updates
    def update_position(self,bounds):
        for i in range(0,num_dimensions):
            self.position_i[i]=self.position_i[i]+self.velocity_i[i]

            # adjust maximum position if necessary
            if self.position_i[i]>bounds[i][1]:
                self.position_i[i]=bounds[i][1]

            # adjust minimum position if neseccary
            if self.position_i[i]<bounds[i][0]:
                self.position_i[i]=bounds[i][0]


def NichePSO(fitness_function,lst_x0,bounds,maxiter,th_,theta_s=False,x_s=False):
    # sub swarm radius
    sub_swarm_rad = 1

    # Params PSO
    w = 0.5     # constant inertia weight (how much to weigh the previous velocity)
    c1 = 1.      # cognative constant
    c2 = 2.      # social constant

    # err_best_g=-1 # best error for group
    # pos_best_g=[] # best position for group
    lst_err_best_g = []
    lst_pos_best_g = []

    # main swarm
    main_swarm = []

    # position initial of particles
    num_particles = len(lst_x0)

    # create particles
    lst_particles=[]
    for i in range(0,num_particles):
        lst_particles.append(Particle(lst_x0[i]))

    # reserved niches
    reserved_niches = []
    reserved_niches_err = []

    # new bounds
    lst_bounds = []

    #
    # Start algorithm NichePSO ...
    #
    for ite in range(maxiter):
        if ite == 0:
            # get first the position of the best particles to help a NichePSO
            main_swarm, lst_err_best_g, lst_pos_best_g = get_fuzzy_niches(lst_particles, fitness_function, th_)

            # update dst between particles
            lst_sub_swarm_rad = []
            for i in range(len(lst_pos_best_g) - 1):
                p0 = lst_pos_best_g[i]
                p1 = lst_pos_best_g[i + 1]
                if len(p0) > 0 and len(p1) > 0:
                    lst_sub_swarm_rad.append(distance.euclidean(p0, p1))
            if len(lst_sub_swarm_rad) > 0:
                sub_swarm_rad = np.min(lst_sub_swarm_rad)
                sub_swarm_rad /= 4.

            # join nearest particles to each niche
            lst_x0, main_swarm = get_particles_nearest(main_swarm, lst_x0, sub_swarm_rad)
            for i in range(0,num_particles):
                lst_particles[i].position_i = list(lst_x0[i])

            # create new bounds
            for i in range(len(lst_pos_best_g) - 1):
                if len(lst_pos_best_g[i]) > 0:
                    x_ = lst_pos_best_g[i][1]
                    lst_bounds.append([(0, 0), (x_ - sub_swarm_rad, x_ + sub_swarm_rad)])
                else:
                    lst_bounds.append([])
            lst_bounds.append(bounds)

            # train PSO algorithm only with cognitive constant
            lst_prts = [lst_particles[k] for k in main_swarm[-1]]
            nss = len(lst_prts)
            lst_prts, _, _ = PSO(nss, bounds, reserved_niches, lst_prts, fitness_function, lst_err_best_g[-1], lst_pos_best_g[-1], th_, w=w, c1=c1, c2=0)
            for k in range(nss):
                lst_particles[main_swarm[-1][k]] = lst_prts[k]

            """
            print "-------"
            print "Zero main swarm"
            print main_swarm
            print "New bounds"
            print lst_bounds
            print "new dst between particles"
            print sub_swarm_rad
            print "first pos best"
            print lst_pos_best_g
            print "first error best"
            print lst_err_best_g
            """

            # get position from particles
            # lst_x0 = get_pos_particles_from_lst(lst_particles)
            # aux_positions_particles = get_particles(lst_x0, np.copy(main_swarm[-1]))
            # main_swarm, lst_err_best_g, lst_pos_best_g = add_new_sub_swarms(main_swarm, aux_positions_particles, lst_err_best_g, lst_pos_best_g)


            # algorithm to create niches
            # main_swarm, lst_err_best_g, lst_pos_best_g = create_sub_swarms(lst_x0)

            # get the  best of the without sub swarm ... in the init ...
            # main_swarm, idx_best_errori, idx_best_positionsi = get_best_of_without_sub_swarm(main_swarm, lst_particles, fitness_function, th_)
            # idx_best_errori, idx_best_positionsi = remove_repeat(idx_best_errori, idx_best_positionsi, step_theta_=0.2)

            # lst_x0_t = np.transpose(lst_x0)
            # plt.scatter(lst_x0_t[0], lst_x0_t[1], s=3.0)
            # plt.show()

        else:
            # Get and return some particles ...
            for i in range(len(main_swarm)):
                sub_swarm = main_swarm[i]
                nss = len(sub_swarm)
                if nss > 1:
                    lst_prts = [lst_particles[k] for k in sub_swarm]
                    lst_prts, lst_err_best_g[i], lst_pos_best_g[i] = PSO(nss, lst_bounds[i], reserved_niches, lst_prts, fitness_function, lst_err_best_g[i], lst_pos_best_g[i], th_, w=w, c1=c1, c2=c2)

                    for k in range(nss):
                        lst_particles[sub_swarm[k]] = lst_prts[k]

            # get position from particles
            lst_x0 = get_pos_particles_from_lst(lst_particles)

        # 1 -> without_sub_swarm to sub_swarm
        # 2 -> sub_swarm to without_sub_swarm
        # main_swarm = check_particle_sub_swarms(main_swarm, lst_x0, opt=2)
        # main_swarm = check_particle_sub_swarms(main_swarm, lst_x0, opt=1)

        #main_swarm = check_merge_sub_swarms(main_swarm, lst_x0)
        #main_swarm, lst_err_best_g, lst_pos_best_g = send_lonely_particles_to_main_swarm(main_swarm, lst_err_best_g, lst_pos_best_g)
        #main_swarm, lst_err_best_g, lst_pos_best_g = clean_sub_swarms(main_swarm, lst_err_best_g, lst_pos_best_g)

    """
    print ("Ite : %s" % (ite))
    print "--------"
    print "main swarm"
    print main_swarm
    print "Positions"
    print lst_x0
    print "Error Best"
    print lst_err_best_g
    print "Pos Best"
    print lst_pos_best_g
    print sub_swarm_rad
    """

    # get the  best of the without sub swarm ... in the final ...
    #main_swarm, idx_best_errorf, idx_best_positionsf = get_best_of_without_sub_swarm(main_swarm, lst_particles, fitness_function, th_)
    #idx_best_errorf, idx_best_positionsf = remove_repeat(idx_best_errorf, idx_best_positionsf, step_theta_=0.2)

    # get reserved niches ...
    reserved_niches_err, reserved_niches = remove_repeat(lst_err_best_g, lst_pos_best_g, step_theta_=sub_swarm_rad)
    reserved_niches_err, reserved_niches = filter_to_minimal(reserved_niches_err, reserved_niches, error)


    """
    print "----"
    print idx_best_errori
    print idx_best_positionsi
    print "----"
    print idx_best_errorf
    print idx_best_positionsf
    print "Reserved niches"
    print reserved_niches_err
    print reserved_niches
    print ""
    """

    # try to plot ...
    #if theta_s is not False and x_s is not False:
    #plot_niche_pso(lst_x0, main_swarm, th_, theta_s, x_s)
    #plot_niche_pso2(False, th_, theta_s, x_s, reserved_niches)

    return reserved_niches_err, reserved_niches

# general function
def function_to_val(x, c):
    #return np.polyval([-1, 0, 1, 0, c, 0], x)
    #return np.polyval([-1, 0, 4, c], x)
    return np.polyval([-1, 1, c], x)

#
# Tools for PSO Algorithm
if __name__ == "__main__":

    # num dimentions
    num_dimensions = 2

    #
    # Tools for function to optimize

    # params
    min_theta1 = -10
    max_theta1 = 10
    min_x = -10
    max_x = 10

    # Warning : update params
    step_theta1 = 0.001     # 1e-3
    step_x = 0.001          # 1e-3
    error_th_m = 0.012             

    # update params
    n_step_r = 0.01
    max_x += step_x
    max_theta1 += step_theta1


    # ranges
    theta2 = np.arange(min_theta1, max_theta1 + n_step_r, n_step_r)  # get a range for theta
    x_s = np.arange(min_x, max_x, step_x)

    reject_particles_error = 1e-12

    #
    # Other params for this algorithm
    n_particles_niche = 260
    th_bisection = 1e-6
    maxiter = 65                 # max number of iter
    dst_max_particles = 0.0008   # max distance between particles for euclidean distance
    error = 1e-3  # -6

    #
    # Call the algorithm
    #
    error_th = .0001
    mtrx_reserved_niches = []  # struct of this matrix [[num roots, value theta, niches, error niches], [], ..]

    # get niches for each
    total_time = 0.
    a = time.time()
    for i in range(len(theta2)):
        th2 = theta2[i]

        # get theta values
        theta_s = [function_to_val(x_, th2) for x_ in x_s]

        # bouns
        first_bounds = [(0, 0), (min_x, max_x)]

        # first run NichePSO Algorithm ...
        lst_pos_particles = create_particles_position(first_bounds, type_bound=2, n_rows=n_particles_niche, n_cols=0.)
        # plot_niche_pso2(lst_pos_particles, th_, theta_s, x_s, reserved_niches)

        # call niche pso algorithm
        reserved_niches_err, reserved_niches = NichePSO(fitness, lst_pos_particles, first_bounds, maxiter, th2, theta_s, x_s)

        # call nitche pso modified algorithm
        # FIXME ...

        # get new niches
        mtrx_reserved_niches.append([len(reserved_niches), th2, reserved_niches, reserved_niches_err])

    b = time.time(); d = b - a; total_time += d; print ("Time niche pso : %s" % (d))
    print ("Num of nitches with Niche PSO : %s " % (len(mtrx_reserved_niches)))
    plt_niche_pso3(mtrx_reserved_niches)
    #print_mtrx(mtrx_reserved_niches)

    # explore new zones with bisection
    a = time.time()
    #mtrx_reserved_niches = check_explore_with_bisection(mtrx_reserved_niches)
    #print_mtrx(mtrx_reserved_niches)
    print ("Num of nitches with bisection method : %s " % (len(mtrx_reserved_niches)))
    b = time.time(); d = b - a; total_time += d; print ("Time bisection method : %s" % (d))

    # create a difurcation diagram
    a = time.time()
    lst_red = []
    lst_blue = []
    for i in range(len(mtrx_reserved_niches) - 1):
        new_route = []
        reserved_niches = []
        reserved_niches_err = []

        nroots1 = mtrx_reserved_niches[i][0]
        nroots2 = mtrx_reserved_niches[i + 1][0]

        nth1 = mtrx_reserved_niches[i][1]
        nth2 = mtrx_reserved_niches[i + 1][1]


        if nroots1 <= nroots2:
            new_route = np.arange(nth2, nth1, - step_theta1)
            reserved_niches = list(mtrx_reserved_niches[i + 1][2])
            reserved_niches_err = list(mtrx_reserved_niches[i + 1][3])

        elif nroots1 > nroots2:
            new_route = np.arange(nth1, nth2 + step_theta1, step_theta1)
            reserved_niches = list(mtrx_reserved_niches[i][2])
            reserved_niches_err = list(mtrx_reserved_niches[i][3])        

        num_reserved_niches = len(reserved_niches)
        history_follow_particles = [[] for _ in range(num_reserved_niches)]
        for nth_ in new_route:
            for k in range(num_reserved_niches):
                rn = reserved_niches[k]
                if len(rn) > 0:
                    # use a umbral to converge only with five points around the center 
                    if len(history_follow_particles[k]) >= 2:  # >= 2
                        m_point = get_ms(history_follow_particles[k], points=2)
                        lst_pos_particles = create_particles_position([0, rn[1]], type_bound=0, n_rows=5., n_cols=0., threshold=error_th_m, m_point=m_point)
                        bounds = [(0, 0), (rn[1] - error_th_m, rn[1] + error_th_m)]
                    else:
                        lst_pos_particles = create_particles_position([0, rn[1]], type_bound=1, n_rows=8., n_cols=0., threshold=error_th)            
                        bounds = [(0, 0), (rn[1] - error_th, rn[1] + error_th)]

                    num_particles = len(lst_pos_particles)
                    lst_particles = []
                    lst_particles_u=[]
                    lst_particles_v=[]
                    for l in range(num_particles):
                        lst_particles.append(Particle(lst_pos_particles[l]))
                        lst_particles_u.append(Particle(lst_pos_particles[l]))
                        lst_particles_v.append(Particle(lst_pos_particles[l]))

                    # call only pso
                    err_best_g = -1
                    pos_best_g = []

                    # call differential pso algorithm
                    #for l in range(30):
                    #    lst_particles, err_best_g, pos_best_g = PSO(num_particles, bounds, [rn], lst_particles, fitness, err_best_g, pos_best_g, nth_, w=0.5, c1=0.8, c2=2)


                    # call diferential evolution
                    err_best_g, pos_best_g = differential_evolution(bounds, lst_particles, lst_particles_v, fitness, nth_, num_dimensions, maxiter=30)

                    reserved_niches[k] = list(pos_best_g)
                    reserved_niches_err[k] = err_best_g

                    # call diferential evolution


                    #print reserved_niches
                    #print reserved_niches_err
                    #print "---------------"

                    # create a bifurcation diagram
                    reserved_niches_aux = abs(reserved_niches[k][1])
                    if reserved_niches_aux < 0.7071 and reserved_niches_aux > 2 * step_theta1:
                        lst_red.append([nth_, reserved_niches[k][1]])

                    elif reserved_niches_aux >= 0.7071:
                        lst_blue.append([nth_, reserved_niches[k][1]])

                    elif reserved_niches_aux >= 0 and reserved_niches_aux <= 2 * step_theta1:
                        if nth_ <= 0:
                            lst_blue.append([nth_, reserved_niches[k][1]])
                        else:
                            lst_red.append([nth_, reserved_niches[k][1]])

                    # update history particles with the last 3 particles ...
                    if len(history_follow_particles[k]) == 2:
                        history_follow_particles[k].pop(0)
                    history_follow_particles[k].append(reserved_niches[k])

                    # update the general list of particles
                    reserved_niches_err, reserved_niches = remove_repeat(reserved_niches_err, reserved_niches, step_theta_=0.0005)
                    reserved_niches_err, reserved_niches = filter_to_minimal(reserved_niches_err, reserved_niches, 0.0045)  # 0.00018

                    num_reserved_niches2 = len(reserved_niches)
                    if num_reserved_niches2 != num_reserved_niches:
                        # delete history of particles not used
                        lst_idx = []
                        for idx in range(num_reserved_niches - 1, -1, -1):
                            if len(history_follow_particles[idx]) > 0:                   
                                last_point_h = history_follow_particles[idx][-1]
                                try:
                                    reserved_niches.index(last_point_h)
                                except:
                                    lst_idx.append(idx)
                        for idx in lst_idx:
                            history_follow_particles.pop(idx)

                        for _ in range(abs(num_reserved_niches2 - num_reserved_niches)):
                            reserved_niches_err.append(-1)
                            reserved_niches.append([])
                            history_follow_particles.append([])

    b = time.time(); d = b - a; total_time += d; print ("Time only pso : %s" % (d))
    print (":> Total time : %s <:" % (total_time))

    lst_blue_t = np.transpose(lst_blue)
    lst_red_t = np.transpose(lst_red)

    plt.scatter(lst_blue_t[0], lst_blue_t[1], s=0.8, color="blue")  # blue
    plt.scatter(lst_red_t[0], lst_red_t[1], s=0.8, color="blue")  # red
    #plt.xlim(-0.4, 0.4); plt.ylim(-1.5, 1.5);
    plt.title("Niche PSO Algorithm", fontsize=16.)
    plt.xlabel('Theta axis'); plt.ylabel('X axis');
    plt.grid()
    plt.show()
