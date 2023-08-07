#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 14:22:13 2023

@author: Cate
"""
import numpy as np
import random
import math

class Nuclei(object):
    

    def __init__(self, decayConstant, length, timestep):
        """
        Initialize a Nuclei with 3 required parameters
        
        Parameters
        ----------
        decayConstant : Double
        length : Int
        timestep : Double


        """
        self.array = np.ones((length, length))
        self.decayConstant = decayConstant
        self.length = length
        self.timestep = timestep

    def radioactiveDecay(self):
        
        """
        
        Returns
        -------
        Estunated half live after simulated radioactive decay

        """

        time = 0
        initialundecayed = self.length * self.length
        undecayed = initialundecayed
        goalDecayed = undecayed/2
        while undecayed > goalDecayed:
            probability =  self.decayConstant * self.timestep
            for i in range(self.length):
                for j in range(self.length):
                    nuclei = self.array[i][j]
                    if nuclei == 1:
                        randomNumber = random.uniform(0,1)
                        if probability > randomNumber:
                            #
                            undecayed -= 1
                            self.array[i][j] = 0
                      
            time += self.timestep
          
        print("Initial number of undecayed nuclei: " + str(initialundecayed))
        print("Final number of undecayed nuclei: " + str(undecayed))
        print("Simulated value of half life: " + str(time))
        calculatedHalfLife = math.log(2, math.e) / self.decayConstant
        print("Actual value of half life: " + str(calculatedHalfLife))
        Grid2D(self.array).printGrid()
        
            
       
class Grid2D(object):
    def __init__(self, arr):
        """
        Grid2D object created with the array given in parameter
        """
        self.grid = arr
        
    def printGrid(self):
        """
        Prints out the 2D grid version of the 2D array
        """
        display = ""
        for i in range(len(self.grid)):
            for j in range(len(self.grid[i])):
                display += str(self.grid[i][j]) + " "
            display += "\n"
        print(display)


def main():
    obj = Nuclei(0.02775, 50, 0.01)
    obj.radioactiveDecay()
    
main()