#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 14:35:44 2023

@author: Cate
"""
import numpy as np
import matplotlib.pyplot as plt

class Mandelbrot(object):
    
    def __init__(self):
        """
        ititialize the mandelbrot set to have the real number range and imaginary number range for the plot

        Returns
        -------
        None.

        """
        
        self.xStart = -2.025
        self.xEnd = 0.6
        self.yStart = -1.125
        self.yEnd = 1.125
        
    
    def draw(self, threshhold = 255):
        """
        function to plot the mandelbrot set. can change threshhold, but automatically set to 255

        Parameters
        ----------
        threshhold : TYPE, integer
            DESCRIPTION. The default is 255. Determines the threshold for iterations when determining if 
            point is in mandelbrot set.

        Returns
        -------
        None. Outputs a colored 2D grid of the mandelbrot set.

        """
        #generate evenly spaced out array of 255 points for x and y of the complex nimber c
        xArray = np.linspace(self.xStart, self.xEnd, threshhold)
        yArray = np.linspace(self.yStart, self.yEnd, threshhold)
        
        #create bounds for the imshow graph
        left = xArray.min()
        right = xArray.max()
        bottom = yArray.min()
        top = yArray.max()
        
        #create an empty array to put the corresponding results from mandelbrot function into
        iterations = np.empty((len(xArray), len(yArray)))

        #iterate through the entire plot and find the value of iterations for each complex number C
        for i in range(len(xArray)):
            for j in range(len(yArray)):
                c = complex(xArray[i], yArray[j])
                #flipped x and y for how numbers represented in matricies
                iterations[j,i] = self.mandelbrot(c)
               
        #create the meshgrid with the real values on X and imaginary values on Y axis
        #xx, yy = np.meshgrid(xArray, yArray)
        
        #show the heatmap according to the number of iterations the mandelbrot set took
        #extent for x and y ticks
        plt.imshow(iterations, extent = (left, right, bottom, top), interpolation = 'none', cmap = 'hot') 
        
        
        
        plt.colorbar()

        plt.show()





    def mandelbrot(self, c):
        """
        Function to determine iterations needed for point c to diverge. Determine if point is in the set. 

`
        Parameters
        ----------
        c : complex number
            point on the plot used to determine iteration/N number.

        Returns
        -------
        iterations : Integer
            The number of iterations it took for c to diverge or reach threshhold.

        """
        #Z always starts at zero as initial value
        z = complex(0,0)
        
        #iterations (N) determine the depth of the color of plotted point
        iterations = 0
        
        #check to see if z has diverged or if iterations have reached threshold. 
        while iterations < 255 and abs(z) < 2: 
            z = z**2 + c
            iterations += 1
            
        #return iterations/N to be used to color the point c on graph
        return iterations

    
        
def main():
    test = Mandelbrot()
    test.draw()
    
main()
                        
                    
                    
              
     
     
        
        