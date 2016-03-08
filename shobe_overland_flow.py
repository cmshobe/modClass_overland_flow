# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 10:08:11 2016

@author: Charlie Shobe

Modeling class: 1-D overland flow with unit test

Models hillslope hydrology under a sinusoidal forcing.

Plots hillslope profile, zoomed-in reach profile, and discahrge at base of hillslope.

"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys

matplotlib.rcParams.update({'font.size': 16})

class Overland(object):
    
    def law_of_wall_unit_test(self, g, depth, slope, z_0): #unit test for LOTW
        print 'Ececuting unit test...'
        test_v = overland_flow.calc_avg_velocity(g, depth, slope, z_0)
        if np.isclose(test_v, 1.434157, 10e-5, 10e-5) == 1:
            print 'Unit test successful'
        else:
            sys.exit('UNIT TEST FAILED')

    def calc_avg_velocity(self, g, h, bed_slope, z_0): #law of the wall
        avg_velocity = np.sqrt(g * h * bed_slope) / .408 * (np.log(h / z_0) - 1)
        return avg_velocity
        
    def plot(self, flow, close_up, discharge, x, z_bedrock, times, discharges, current_time, it, z_total):
        flow.clear()
        close_up.clear()
        discharge.clear()
        flow.plot(x, z_bedrock, color='k', linewidth = 2, label='Bedrock')
        flow.plot(x, z_total, color='b', linewidth = 2, label='Water')
        flow.set_xlim(0, 100)
        flow.set_ylim(min(z_bedrock)+5, max(z_bedrock)+5)
        flow.plot((80, 80), (min(z_bedrock)+5, max(z_bedrock)+5), linestyle='--', linewidth=3, color='k')
        flow.plot((85, 85), (min(z_bedrock)+5, max(z_bedrock)+5), linestyle='--', linewidth=3, color='k')
        flow.set_xlabel('Distance [m]')
        flow.set_ylabel('Elevation [m]')
        flow.text(5, (min(z_bedrock)+max(z_bedrock))/2, 'Time [s]: %.1f' % current_time)
        close_up.plot(x, z_bedrock, color='k', linewidth = 2, label='Bedrock')
        close_up.plot(x, z_total, color='b', linewidth = 2, label='Water')
        close_up.set_xlim(80, 85)        
        close_up.set_ylim(0, 10)  
        close_up.set_ylabel('Elevation [m]')
        close_up.set_xlabel('Distance [m]')
        discharge.plot(times[0:it], discharges[0:it], marker='+', color = 'r', linewidth = 0)
        discharge.set_xlabel('Time (s)')
        discharge.set_ylabel('Discharge [m2/s]')
        discharge.set_ylim(0, 0.5)
        discharge.set_xlim(0, max(times))
        plt.subplots_adjust(hspace=0.5)
        plt.pause(0.0001)
    
    def main(self, bed_slope, z_0):
        g = 9.81 #acceleration due to gravity
        rain_mean = 0.0005 #m/s
        infiltration = 0.00001 #m/s
        
        #spatial domain
        x_min = 0 #m
        dx = 0.2 #m
        x_max = 100 #m
        x = np.arange(x_min, x_max + dx, dx)
        
        z_bedrock = 100 - bed_slope * x #this line for inclined plane of slope bed_slope
        p = 10
        k_shape = 60
        h_param = 0
        z_bedrock = ((-(x-h_param)**2)/(12*p))+k_shape + 2
        
        #time domain
        t_min = 0 #s
        dt = .1 #s
        t_max = 10000 #s
        times = np.arange(t_min, t_max + dt, dt)
        t_plot = 100 #plot every __ s
        
        #rainfall distribution is a sinusoid
        rain_distribution = rain_mean + (0.002*np.sin(2*np.pi*times/5000))

        #depth array
        h = 0. * np.ones((len(x))) #for beginning of loop
        discharges = np.zeros((len(times)))
        
        #plotting
        flow_fig = plt.figure(figsize=(14,10)) #instantiate figure
        flow = plt.subplot(311)
        close_up = plt.subplot(312)
        discharge = plt.subplot(313)
        plt.ion()
        plt.show()
    
        ###############FORWARD EULER SOLUTION
        it = 0
        for time in range(len(times)):
            it += 1
            current_time = times[time]
            rain = rain_distribution[it - 1]
            if np.any(h <= 0): #stops ln(0) error in law of the wall
                avg_velocity = 0
            else:
                avg_velocity = overland_flow.calc_avg_velocity(g, h, bed_slope, z_0)
            q = avg_velocity * h
            q = np.insert(q, 0, 0)
            discharges[time] = q[-1]    
            dq_dx = np.diff(q) / dx
            h += (-dq_dx + rain - infiltration) * dt
            h = h.clip(min = 0)
            z_total = z_bedrock + h
            if current_time % t_plot == 0: #plot stuff
                overland_flow.plot(flow, close_up, discharge, x, z_bedrock, times, discharges, current_time, it, z_total)
            else:
                pass
            
        return discharges[-1]
    
if __name__ == "__main__":
    overland_flow = Overland()
    overland_flow.law_of_wall_unit_test(9.81, 1, 0.001, 0.001)
    overland_flow.main(0.001, 0.001)