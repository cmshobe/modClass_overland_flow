# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 09:20:04 2016

@author: Charlie
"""


#from scipy.sparse import spdiags
#from scipy.sparse.linalg import spsolve

#def make_advection_matrices(z, r):
#    """Return matrices A and M for advection equations"""
#    ones = np.ones(len(z))
#    A = spdiags( [-beta*r, ones, beta*r], (-1,0,1), len(z), len(z) )
#    M = spdiags( [(1-beta) * r, ones, -(1-beta) * r], (-1,0,1), len(z), len(z) )
#    return A.tocsr(), M.tocsr()

    
    ###########BRAUN/WILLETT SOLUTION (DOESN'T WORK)
    #for time in range(len(times)):
    #    current_time = times[time]
    #    for cell in range(len(x) - 2, -1, -1):
    #        h[cell] = (h[cell] + ((a * 1.5 * np.power(h[cell], 0.5) * h[cell + 1] * dt) / dx)) / (1 + ((a * 1.5 * np.power(h[cell], 0.5) * dt) / dx))
    #    if current_time % t_plot == 0: #plot stuff
    #        flow.clear()
    #        #flow.plot(x / 1000, zb, color='k', linewidth = 2, label='Bedrock')
    #        flow.plot(x, h, color='b', label='Water')
    #        flow.set_xlim(0, x_max)
    #        flow.set_ylim(0, 1)
    #        plt.xlabel('Distance [m]')
    #        flow.set_ylabel('Elevation [m]')
    #        flow.text(5, .5, 'Time [s]: %.1f' % current_time)
    #        plt.pause(0.1)
    #        #glacier_fig.savefig('glacier'+str(current_time)+'.png')
    #    else:
    #        pass
    #    
    ##########CRANK-NICHOLSON SOLUTION (BROKEN AS HELL)
    # Set up basic constants
    #beta = 0.5
    #J = 200 # total number of mesh points
    #z = np.linspace(-10,10,J) # vertices
    #dz = abs(z[1]-z[0]) # space step
    #dt = 0.001    # time step
    #v = 2 * np.ones(len(z)) # velocity field (constant)
    #r = v / 2 * dt / dx
    
    # Initial conditions (peak function)
    #gaussian = lambda x, height, position, hwhm: height * np.exp(-np.log(2) * ((x - position)/hwhm)**2)
    #h_init = gaussian(x, 1, -3, 2)
    #h[:] = 0.05
    #h[0:50] = np.arange(0.05, 0.1, 0.001)
    #h[50:100] = np.arange(0.1, 0.05, -0.001)
    ##h[:] = h_init
    ##h[0:50] = 0.5
    #for time in range(len(times)):
    #    current_time = times[time]
    #    v[:] = np.sqrt(g * bed_slope / chezy_roughness_coefficient) * 1.5 * np.power(h, 0.5)
    #    r = v / 2 * dt / dx
    #    A, M = make_advection_matrices(h, r)    
    #    h = spsolve(A, M * h)
    #    if current_time % t_plot == 0: #plot stuff
    #        flow.clear()
    #        #flow.plot(x / 1000, zb, color='k', linewidth = 2, label='Bedrock')
    #        flow.plot(x, h, color='b', label='Water')
    #        flow.set_xlim(0, 100)
    #        flow.set_ylim(0, 0.5)
    #        plt.xlabel('Distance [m]')
    #        flow.set_ylabel('Elevation [m]')
    #        flow.text(5, .5, 'Time [s]: %.1f' % current_time)
    #        plt.pause(0.001)
    #    else:
    #        pass
    