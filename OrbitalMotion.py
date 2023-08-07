#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 14:23:58 2023

@author: Cate
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math

class Planet(object):
    def __init__(self, mass, position, initialvelocity):
        self.mass = mass
        self.position = position
        self.initialvelocity = initialvelocity
        
    def getMass(self):
        return self.mass
    
    def getPosition(self):
        #array
        return self.position
    
    def setPosition(self, newposition):
        self.position = newposition
        
    def getVelocity(self):
        return self.initialvelocity
        
    
        
class OrbitalMotion(object):
    def __init__(self, planet_one, planet_two, num_timesteps, length_timestep):
        self.planet_one = planet_one
        self.planet_two = planet_two
        self.num_timesteps = num_timesteps
        self.length_timestep = length_timestep
        
        self.totalx_p1 = []
        self.totaly_p1 = []
        self.totalx_p2 = []
        self.totaly_p2 = []

    
    def simulate(self):
        r1 = np.array(self.planet_one.getPosition())
        r2 = np.array(self.planet_two.getPosition())
        
        velocity_one = self.planet_one.getVelocity()
        velocity_two = self.planet_two.getVelocity()
        self.totalx_p1.append(r1[0])
        self.totaly_p1.append(r1[1])
        
        for x in range(self.num_timesteps):
            print(r1)
            print(r2)
            r12 = r2 - r1
            r21 = r1 - r2
            
            acceleration_one = -6.6743e-11 * (self.planet_two.getMass() / math.pow(np.linalg.norm(r21),3)) * r21
            acceleration_two = -6.6743e-11 * (self.planet_one.getMass() / math.pow(np.linalg.norm(r12),3)) * r12
            
            velocity_one = velocity_one + (acceleration_one * self.length_timestep)
            velocity_two = velocity_two + (acceleration_two * self.length_timestep)
            
            r1 = r1 + velocity_one * self.length_timestep
            r2 = r2 + velocity_two * self.length_timestep
           
            
            kineticenergy = np.linalg.norm(0.5 * self.planet_one.getMass() * velocity_one**2 + 0.5 * self.planet_two.getMass() * velocity_two**2) 
            print("KINETIC ENERGY: " + str(kineticenergy))
            
            self.totalx_p2.append(r2[0])
            self.totaly_p2.append(r2[1])
            self.totalx_p1.append(r1[0])
            self.totaly_p1.append(r1[1])

          
            
    def getx1(self):
        return self.totalx_p1
    def gety1(self):
        return self.totaly_p1
    def getx2(self):
        return self.totalx_p2
    def gety2(self):
        return self.totaly_p2
    
      

# class Animation(object):

#     def __init__(self, xpos, ypos):
#         # create arrays of x and y coordinates for plot
#         self.xpos = xpos
#         self.ypos = ypos

#     def animate(self, i):
#         self.patch.center = (self.xpos[i], self.ypos[i])
#         return self.patch,

#     def run(self):
#         # create plot elements
#         fig = plt.figure()
#         ax = plt.axes()

#         # create circle of radius 0.1 centred at initial position and add to axes
#         self.patch = plt.Circle((self.xpos[0], self.ypos[0]), 1000, color = 'g', animated = True)
#         ax.add_patch(self.patch)

#         # set up the axes
#         ax.axis('scaled')
#         ax.set_xlim(np.min(self.xpos), np.max(self.xpos))
#         ax.set_ylim(np.min(self.ypos), np.max(self.ypos))
#         ax.set_xlabel('x')
#         ax.set_ylabel('y')

      
#         # create the animator
#         self.anim = FuncAnimation(fig, self.animate, frames = len(self.xpos), repeat = False, interval = 1, blit = True)

#         # show the plot
#         plt.show()
        
class PlanetAnimation(object):

    def __init__(self, orbital):
        # set initial and final x coordinates of circle
       
        self.xpos1 = orbital.getx1()
        self.xpos2 = orbital.getx2()
        
        self.ypos1 = orbital.gety1()
        self.ypos2 = orbital.gety2()
        

        # set up simulation parameters
        # self.niter = 500
        # self.xincr =  (np.maxself.xmax - self.xpos)/self.niter
        
    def animate(self, i):
        # update the position of the circles

        self.patch.center = (self.xpos1[i], self.ypos1[i])
        return self.patch,
        # self.patches[0].center = (self.xpos1[i], self.ypos1[i])
        # self.patches[1].center = (self.xpos2[i], self.ypos2[i])
     
        # return self.patches,


    def run(self):
        
        fig = plt.figure()
        ax = plt.axes()
        maxOrb = math.sqrt(np.dot(1.2,1.2))  

        # create circle of radius 0.1 centred at initial position and add to axes
        self.patch = plt.Circle((self.xpos1[0], self.ypos1[0]), 0.05*maxOrb, color = 'g', animated = True)
        # self.patch = plt.Circle((self.xpos1[0], self.ypos1[0]), 10000, color = 'g', animated = True)

        ax.add_patch(self.patch)

        b = 1.2
        lim = maxOrb*b

        ax.axis('scaled')
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)

        


        # set up the axes
        # ax.axis('scaled')
        # # ax.set_xlim(min(np.min(self.xpos2), np.min(self.xpos1)), max(np.max(self.xpos2), np.max(self.xpos1)))
        # ax.set_xlim(min(np.min(self.ypos2), np.min(self.ypos1)), max(np.max(self.ypos2), np.max(self.ypos1)))
        # ax.set_ylim(min(np.min(self.ypos2), np.min(self.ypos1)), max(np.max(self.ypos2), np.max(self.ypos1)))
        
        
        ax.set_xlabel('x (rads)')
        ax.set_ylabel('sin(x)')

        # create the animator
        self.anim = FuncAnimation(fig, self.animate, frames = len(self.xpos1), repeat = False, interval = 1, blit = True)

        # show the plot
        plt.show()
        
        # # create plot elements
        # fig = plt.figure()
        # ax = plt.axes()

        # # create list for circles
        # self.patches = []

        # # create circles of radius 0.1 centred at initial position and add to list
        # self.patches.append(plt.Circle((self.xpos1[0], self.ypos1[0]), 100000, color = 'g', animated = True))
        # self.patches.append(plt.Circle((self.xpos2[0],self.ypos2[0]), 10, color = 'b', animated = True))
        # # add circles to axes
        # for i in range(0, len(self.patches)):
        #     ax.add_patch(self.patches[i])

        # # set up the axes
        # # ax.axis('scaled')
        # ax.set_xlim(min(np.min(self.xpos2), np.min(self.xpos1)), max(np.max(self.xpos2), np.max(self.xpos1)))
        # ax.set_ylim(min(np.min(self.ypos2), np.min(self.ypos1)), max(np.max(self.ypos2), np.max(self.ypos1)))
        # # ax.set_xlim(-100, 9777300)
        # # ax.set_ylim(-10,49106)
        # ax.set_xlabel('x (rads)')
        # ax.set_ylabel('sin(x)')

        # # create the animator
        # self.anim = FuncAnimation(fig, self.animate, frames = len(self.xpos2), repeat = False, interval = 100, blit = True)

        # # show the plot
        # plt.show()

        

def main():
    mars = Planet(6.4185e23, [0,0], [0,0])
   
    velocity_two = np.sqrt(np.abs(-6.6743e-11 * (mars.getMass() / 9.3773e6)))
    phobos = Planet(1.06e16, [9.3773e6, 0], [0,velocity_two])
    motion = OrbitalMotion(phobos, mars, 200, 0.1)
    motion.simulate()
    
    animation = PlanetAnimation(motion)
    animation.run()

main()
        