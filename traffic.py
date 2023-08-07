#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 13:52:11 2023

@author: Cate
"""
import numpy as np
import random
import matplotlib.pyplot as plt


class Traffic(object):
    
    def __init__(self, N):
        """
        

        Parameters
        ----------
        N : Integer. Number of cells/spots on the road

        Returns
        -------
        None.

        """
        self.N = N

   
    
    #INITIALIZE THE ROAD
    def simulate(self, timesteps, density):
        """
        

        Parameters
        ----------
        timesteps: the number of timesteps that the program will cycle through
        
        density: the density of the cars on the road. used with N to calculate number of cars.
        
        Returns
        -------
        None. 
        Prints out the average speed for each time step, the steady state speed for the density.
        Plots the movmement of cars as a function of time.

        """
        numcars = 0 #none before initializing the road
        #density: proportion of cells with cars - initialize with random poisitions
        road = np.zeros(self.N, dtype=int)
        for i in range(self.N):
            rand = random.uniform(0,1)
            if rand <= density:
                road[i] = 1
                numcars += 1 #add car to total cars initialized
        
        #variables for finding the avg speed and steady state speed
        steadystate = 0
        rollingavg = 0
        avgspeed = 0
        
        
        totalmovement = [road]
    
        #empty new road to move cars into
        newroad = np.zeros(self.N, dtype=int)
        
        #loop through all of the time steps
        for t in range(timesteps):
            carsmoved = 0
            
            #at each time step, check if the car can move or should stay in its place
            for j in range(self.N):
                if road[j] == 1: #was 1
                   
                    if road[(j+1) % self.N] == 1:
                        newroad[j] = 1
                    else: 
                        #set this to meaning that a car has moved
                        carsmoved += 1
                        newroad[j] = 0
                else: 
                    if road[(j-1) % self.N] == 1:
                        newroad[j] = 1
                    else: 
                        newroad[j] = 0
                        
            #append this array of car positions to entire history of movement for graph of car movmement
            totalmovement.append(road)
            
            #set the new car positions as the current road to check
            road = np.copy(newroad)
            
            #check numcars value to avoid division by zero if road is empty
            if numcars != 0:
                avgspeed = carsmoved/numcars
            else:
                avgspeed = 1
                
            #once the avgspeed stops changing, set the steadystate speed as that speed. 
            if avgspeed == rollingavg:
                steadystate = avgspeed
            else:
                rollingavg = avgspeed
            print("average speed for timestep: " + str(avgspeed))
        
        
        #plot the car movement. Black square is no car. White square means there is a car.
        plt.imshow(totalmovement, cmap='hot', interpolation='nearest')
        plt.xlabel("cars")
        plt.ylabel("timesteps")
        plt.show()
        
        #once the speed stops changing, this is the stead speed. Print this out along with associated density.
        print("steady state speed for density of " + str(density) + " is " + str(np.round(steadystate, 3)))
        return np.round(steadystate, 3)
  
        
        
            
            
        
def main():
    test = Traffic(10)
    test.simulate(10, 0.9)
    
    #create plot of variaty of densities
    # densities = np.linspace(0,1,10) #array of evenly spaced densities between 0 and 1
    # averageSpeeds = [] 
    # for i in range(len(densities)):
    #     averageSpeeds.append(test.simulate(20, densities[i])) #find the steady state average speeds
    # plt.plot(densities, averageSpeeds) #plot density versus steady state average speed 
    # plt.xlabel("densities")
    # plt.ylabel("average speeds")
    # plt.show()
            
    
main()
        
                    
                