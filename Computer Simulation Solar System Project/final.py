#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Catie Corley
Spring 2023

Computer Simulation Final Project
"""


# matplotlib.use('MacOSX')

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numpy.linalg import norm


'''
Class to run the orbital simulation

'''


class Solar(object):

    def __init__(self):
        self.euler = False
        inputdata = []
        "D:\\DEsk\welcome.txt"
        filename = "parameters-solar.txt"
        filein = open(filename, "r")
        for line in filein.readlines():
            if (not line.startswith("#")):
                inputdata.append(line)
        filein.close()

        # simulation parameters
        self.niter = int(inputdata[0])
        self.dt = float(inputdata[1])
        self.G = float(inputdata[2])

        # list for mars and moons
        self.bodies = []
        earth_radius = 0
        # rest of input data is mars and moon data in four line 'chunks'
        # first entry must be mars
        for i in range(3, len(inputdata)-3, 4):
            name = inputdata[i]
            mass = float(inputdata[i+1])
            orbit = float(inputdata[i+2])
            colour = inputdata[i+3]

            #Get earth orbit for satellite
            if name.strip() == 'earth':
                earth_radius = orbit
                
            self.bodies.append(Planet(name, mass, orbit, colour))
        
        #Create a "satellite" body and add to list of bodies
        satellite = Planet("satellite", 8.89e-20, earth_radius + 0.001, 'brown' )
        self.bodies.append(satellite)
        
        # set initial positions and velocities relative to sun
        # sum must be first element in bodies list!
        for i in range(0, len(self.bodies)):
            self.bodies[i].initialise(self.G, self.bodies[0])
        
        
        #Initialize arrays for energy vs. time graph
        self.totalEnergyArray = []
        self.timestamps = []
        
        #initialize array to write total energy out to file
        self.outputText = []
        
        #Initialize to keep track of doomsday
        self.timeSinceDoomsday = 0
        
        #Variables to keep track of Satellite's distance to Mars
        self.distToMars = 10
        self.satelliteLandingTime = 0
        #Variables to keep track of Satellite's distance back to Earth
        self.distToEarth = 10
        self.satelliteReturnTime = 0


    def init(self):
        # initialiser for animator
        return self.patches


    def animate(self, i):
        '''
        Animate function called every iteration by animation in run() function
        Updates the position and velocity of the bodies
        Checks the planet's alignment, system's total energy, and position of launched satellite
        '''

        # keep track of time in earth years
        time = (i+1)*self.dt
       
        if self.euler:
            # update positions for euler
            for j in range(0, len(self.bodies)):
                self.bodies[j].updatePosEuler(self.G, self.dt)
                self.patches[j].center = self.bodies[j].r
        else:
            
            # update positions
            for j in range(0, len(self.bodies)):
                self.bodies[j].updatePos(self.G, self.dt)
                self.patches[j].center = self.bodies[j].r
            
        if self.euler:
            # then update velocities
            for j in range(0, len(self.bodies)):
                for k in range(0, len(self.bodies)):
                    if (j != k):
                        self.bodies[j].updateVelEuler(self.G, self.dt, self.bodies[k])
        else:
            # then update velocities
            for j in range(0, len(self.bodies)):
                for k in range(0, len(self.bodies)):
                    if (j != k):
                        self.bodies[j].updateVel(self.G, self.dt, self.bodies[k])

        # check year and print year if new year for any planet
        for j in range(0, len(self.bodies)):
            if (self.bodies[j].newYear()):
               print(self.bodies[j].name.strip() + " " + str(self.bodies[j].year) + " years = " + str(time) + " earth years")
               # in new year is earth year, also print total energy
               
               #Calculate and print the orbital period every new year for each of the planets
               orbitalPeriod = str((time/self.bodies[j].year) * 365)
               print("Orbital period for " + self.bodies[j].name + " = " + str(orbitalPeriod) + " days.")
               
               if (self.bodies[j].name.strip() == 'earth'):
                   # need to convert from earth masses AU^2 yr^-2 to kg m^2 s-2 (J)
                   # 1 earth mass = 5.97219e24 kg
                   # 1 AU = 1.496e+11 m
                   c =(5.97219e+24*1.496e+11*1.496e+11)/(3.154e+7*3.154e+7)
                   energy = self.energy()*c
                   print('Time = ' + str(time) + ' earth years. Total energy = ' + '{:.3e}'.format(energy) + ' J')
       
        #Check if doomsday occured
        if self.doomsdayCheck():
            #if doomsday is true and planets have aligned, set time since last doomsday to 0
            self.timeSinceDoomsday = 0
        else:
            #no doomsday has occured to add timestep to time since doomsday
            self.timeSinceDoomsday += self.dt
            
         
        #Use below print statement to track doomsdays
        #print("Time since last doomsday: " + str(self.timeSinceDoomsday))
        
        #Calculate the total energy of the system.
        self.calcTotalEnergy(i)  
        #Plot the total energy of the system versus time so far at exactly 1000 iterations
        if i == 3000:
            self.graphTotalEnergy()
            
        #Track the satellite with below function calls
        # self.trackSatelliteToMars(time)
        # self.trackSatelliteReturn(time)
        
        
        return self.patches

    def trackSatelliteToMars(self, currentTime):
         '''
        Function to track the Satellite with respect to mars
        Keeps track of the closest distance Satellite gets to Mars
        Keeps track of the time that Satellite passes over Mars
        '''
        #Find both Mars and the Satellite in bodies array
         for i in range(len(self.bodies)):
             if (self.bodies[i].name.strip() == 'mars'):
                 mars = self.bodies[i]
             elif (self.bodies[i].name.strip() == 'satellite'):
                 sat = self.bodies[i]
         #check if Satellite is passing over mars
         currentDistancefromMars = self.checkDistanceBetweenPlanets(mars, sat) 
         if currentDistancefromMars < 0.11:
             print("Satellite passing over Mars!")
    
        #Keep track of the closest the satellite ever gets to mars and when this occurs
         if currentDistancefromMars < self.distToMars:
             self.distToMars = currentDistancefromMars
             self.satelliteLandingTime = currentTime
             #Reset the min distance to earth to now track the satellite for landing back on earth
             self.distToEarth = 10 
         
         #print out the minimum distance between Satellite and Mars
         print("Min distance btwn Satellite and Mars: " + str(self.distToMars) + " Occured at: " + str(self.satelliteLandingTime))
            
         
            
    def trackSatelliteReturn(self, currentTime):
         '''
         Function to track the Satellite with respect to Earth
         Keeps track of the closest distance Satellite gets to Earth after passing over Mars
         Keeps track of the time that Satellite returns back to Earth
         '''
        #Find both Earth and the Satellite in bodies array
         for i in range(len(self.bodies)):
             if (self.bodies[i].name.strip() == 'earth'):
                 earth = self.bodies[i]
             elif (self.bodies[i].name.strip() == 'satellite'):
                 sat = self.bodies[i]
         
        #check if Satellite is passing over Earth
         currentDistancefromEarth = self.checkDistanceBetweenPlanets(earth, sat) 
         if currentDistancefromEarth < 0.11:
             print("Satellite returning to Earth!")
    
        
        #Keep track of the closest the satellite ever gets to Earth after passing over Mars and when this occurs
         if currentDistancefromEarth < self.distToEarth:
             self.distToEarth = currentDistancefromEarth
             self.satelliteReturnTime = currentTime
         
         #print out the minimum distance between Satellite and Earth
         print("Min distance btwn Satellite and Earth: " + str(self.distToEarth) + " Occured at: " + str(self.satelliteReturnTime))
            
         
            
         
    def energy(self):
        '''
        File to calculate and return total energy of system
        '''
        ke = 0.0
        pe = 0.0
        for j in range(0, len(self.bodies)):
            ke += self.bodies[j].kineticEnergy()
            for k in range(0, len(self.bodies)):
                if (k != j):
                    r = norm(self.bodies[k].r - self.bodies[j].r)
                    pe -= self.G*self.bodies[j].m*self.bodies[k].m / r
        # divide pe by two to avoid double countin
        pe = pe / 2
        totEnergy = ke + pe
        return totEnergy




    def calcTotalEnergy(self, i):
        '''
        Function to calculate the total energy of the system (kinetic + gravitational potential energy)
        Adds current energy to total energy array
        Writes current system energy to file
        '''
        ke = 0.0
        pe = 0.0
        for j in range(0, len(self.bodies)):
            ke += self.bodies[j].kineticEnergy()
            for k in range(0, len(self.bodies)):
                if (k != j):
                    r = norm(self.bodies[k].r - self.bodies[j].r)
                    pe -= self.G*self.bodies[j].m*self.bodies[k].m / r
        # divide pe by two to avoid double countin
        pe = pe / 2
        totEnergy = ke + pe
        
        #Every 100 iterations, print out the total energy of the system to the terminal
        if i % 100 == 0:  
            print('Time = ' + str(i) + ' iterations. Total energy = ' + '{:.3e}'.format(totEnergy)) 
        
        #Add the time and total energy to arrays used for plotting
        self.totalEnergyArray.append(totEnergy)
        self.timestamps.append(i)   
        
        #Every 50 iterations, the total energy of the system is written to the file "SystemsTotalEnergy.txt"
        if i % 50 == 0:
            self.outputText.append('Time = ' + str(i) + ' iterations. Total energy = ' + '{:.3e}'.format(totEnergy) + '\n')
            file_object = open('SystemsTotalEnergy.txt', 'w')
            file_object.writelines(self.outputText)
            file_object.close()
        
        
        
        
    def graphTotalEnergy(self):
        '''
        Function to graph the total energy of the system versus time
        ''' 
        #create new figure to display graph
        fig = plt.figure()
        sub = fig.add_subplot()
        #plot time vs total energy arrays
        sub.plot(self.timestamps, self.totalEnergyArray)
        #Title based on iteration method
        if self.euler:
            sub.title.set_text('Total Energy vs. Time (Euler)')
        else:
            sub.title.set_text('Total Energy vs. Time (Beeman)')
        #Scale axes appropriately
        sub.set_ylim(np.min(self.totalEnergyArray) - 0.03, np.max(self.totalEnergyArray) + 0.03)
        sub.set_xlabel('Iterations')
        sub.set_ylabel('Total Kinetic Energy')
        fig.show()
        
        
        
        
    def doomsdayCheck(self):
        '''
        Function to check if doomsday has occured.
        Will print out "Doomsday" to terminal when all 5 innermost planets are within 5 
        degrees of the mean angle of the planets (AKA Doomsday)
        '''
        #Adjust all positions so that the sun is at (0,0) 
        adjustedPositions = []
        sunPosition = self.bodies[0].r
        for i in range(1, len(self.bodies )- 1):
            adjustedPositions.append(self.bodies[i].r - sunPosition)
        
        #Find the angle of each planet in degrees and add to the array "angles"
        angles = []
        for i in range(len(adjustedPositions)):
            #Use arcTan to find the angle and convert to degrees
            angle = math.atan(adjustedPositions[i][1] / adjustedPositions[i][0])
            angles.append(math.degrees(angle))
        
        #Find the mean angle of all the angles added to array
        meanAngle = np.mean(angles)

        #Check if any of the planets are greater than 5 degrees away from the mean angle
        #If any planet is outside of 5 degrees of the mean angle, return false (Doomsday has not occured)
        for i in range(len(angles)):
            if np.abs(angles[i] - meanAngle) > 5:
                return False
            
        #Print out doomsday to terminal to signal planetary alignment
        print("DOOMSDAY!!!")
        
        #All planets within 5 degrees of mean angle. Return True to alert that doomsday has occured
        return True                      
    
    
    
    

    def run(self):
        '''
        Function to run the simulation with Beeman method
        '''
        
        # set up the plot components        
        fig = plt.figure()
        ax = plt.axes()

        # create an array for patches (planet and moons)
        self.patches = []

        # get orbital radius of outermost moon to set size of orbiting bodies and of plot
        # hacky - should really check to see which moon is outermost
        
        #changed to -2 to exclude satellite 
        maxOrb = math.sqrt(np.dot(self.bodies[-2].r, self.bodies[-2].r))

        # add the planet and moons to the Axes and patches
        for i in range(0, len(self.bodies)):
            if (i == 0):
                self.patches.append(ax.add_patch(plt.Circle(self.bodies[i].r, 0.05*maxOrb, color = self.bodies[i].c, animated = True)))
            else:
                self.patches.append(ax.add_patch(plt.Circle(self.bodies[i].r, 0.02*maxOrb, color = self.bodies[i].c, animated = True)))
        
        # set up the axes
        # scale axes so circle looks like a circle and set limits with border b for prettier plot
        b = 1.2
        lim = maxOrb*b
        print(lim)
        print(self.niter)
        ax.axis('scaled')
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
       
        ax.title.set_text('Solar System Simulation (Beeman)')


        self.anim = FuncAnimation(fig, self.animate, init_func = self.init, frames = self.niter, repeat = False, interval = 1, blit= True)

        plt.show()
        
        
        
        
    def runEuler(self):
        '''
        Function to run the simulation with Euler method
        '''
        self.euler = True
        # set up the plot components        
        fig = plt.figure()
        ax = plt.axes()

        # create an array for patches (planet and moons)
        self.patches = []

        # get orbital radius of outermost moon to set size of orbiting bodies and of plot
        # hacky - should really check to see which moon is outermost
        
        #changed to -2 to exclude satellite 
        maxOrb = math.sqrt(np.dot(self.bodies[-2].r, self.bodies[-2].r))

        # add the planet and moons to the Axes and patches
        for i in range(0, len(self.bodies)):
            if (i == 0):
                self.patches.append(ax.add_patch(plt.Circle(self.bodies[i].r, 0.05*maxOrb, color = self.bodies[i].c, animated = True)))
            else:
                self.patches.append(ax.add_patch(plt.Circle(self.bodies[i].r, 0.02*maxOrb, color = self.bodies[i].c, animated = True)))
        
        # set up the axes
        # scale axes so circle looks like a circle and set limits with border b for prettier plot
        b = 1.2
        lim = maxOrb*b
        print(lim)
        print(self.niter)
        ax.axis('scaled')
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.title.set_text('Solar System Simulation (Euler Method)')
    
        self.anim = FuncAnimation(fig, self.animate, init_func = self.init, frames = self.niter, repeat = False, interval = 1, blit= True)

        plt.show()

    def checkDistanceBetweenPlanets(self, planet_one, planet_two):
        '''
        Function to calculate the distance between two given planets
            Used for tracking satellite
        '''
        distance = np.sqrt((planet_one.r[0] - planet_two.r[0])**2 + (planet_one.r[1] - planet_two.r[1])**2)
        return distance







'''
Planet class

'''


class Planet(object):
 
    def __init__(self, name, mass, orbit, colour):
        self.name = name
        # mass in kg
        self.m = mass
        # orbital radius in m
        self.orbit = orbit
        # colour - need to strip trailing line return!
        self.c = colour.strip()
        # set year to zero
        self.year = 0

    def initialise(self, G, p):
        # inital position, initial coords = (orbit radius, 0)
        self.r = np.array([self.orbit, 0])
        # inital velocity, tangential to position
        # speed = sqrt(G*marsmass/r)
        if (self.orbit == 0.0):
            self.v = np.array([0, 0])
        elif self.name.strip() == 'satellite':
            #launch the satellite at speed so that it drives by Mars and doesn't perfectly orbit like the planets
            self.v = np.array([2.7, 5.9])
        else:
            vel = math.sqrt(G*p.m/self.orbit)
            self.v = np.array([0, vel])
        # intial accelatation, using gravitational force law
        if (self.orbit == 0.0):
            self.a = np.array([0, 0])
        else:
            self.a = self.updateAcc(G, p)
        # set acc_old = acc to start Beeman
        self.a_old = self.a
        
    
    # Integrating with Direct Euler
    def updatePosEuler(self, G, dt):
        self.r_old = self.r
        self.r = self.r + self.v*dt
        
    def updateVelEuler(self, G, dt, p):
        a_new = self.updateAcc(G, p)
        self.v = self.v + a_new*dt
        
        self.a_old = self.a
        self.a = a_new


    # Integrating with Beeman Method
    def updatePos(self, G, dt):
        # keep old position to check for year
        self.r_old = self.r
        
        # update position first: Beeman
        self.r = self.r + self.v*dt + (4*self.a - self.a_old)*dt*dt/6.0
        
    def updateVel(self, G, dt, p):
        # update velocity second: Beeman
        a_new = self.updateAcc(G, p)
        self.v = self.v + + (2*a_new + 5*self.a - self.a_old)*dt/6.0
        # now update acc ready for next iteration
        self.a_old = self.a
        self.a = a_new

    def updateAcc(self, G, p):
        # update acc (gravitational force law)
        pos = self.r - p.r
        a = -G*p.m*pos/math.pow(norm(pos),3)
        return a


    def newYear(self):
        # update the year when the planet passes the +x axis
        if (self.r_old[1] < 0.0 and self.r[1] >= 0.0):
            self.year +=1
            return True
        else:
            return False


    # determine kinetic energy
    def kineticEnergy(self):
        # ke in J
        ke = (np.dot(self.v, self.v))*self.m/2
        return ke


'''
Run Simulation and Experiments
'''
def main():
    
    #Run experiments using Beeman Method
    s = Solar()
    s.runEuler()
    
    #Run experiments using Euler Method
    r = Solar()
    r.run()
    
main()