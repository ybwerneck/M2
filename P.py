# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 16:58:20 2023

@author: yanbw
"""

import numpy as np
import matplotlib.pyplot as plt

def poission(dl,L):
    # Define the problem parameters
    Lx =L # Length of the domain in x-direction
    Ly = L  # Length of the domain in y-direction
    dx=dy=dl
    Nx = Ny= int(L/dl)
    # Define the boundary conditions
    u_top = 0.0
    u_bottom = 0.0
    u_left = 0.0
    u_right = 0.0
    
    # Create the grid
    x = np.linspace(0, Lx, Nx)
    y = np.linspace(0, Ly, Ny)
    X, Y = np.meshgrid(x, y)
    
    # Initialize the u matrix
    u = np.zeros((Ny, Nx))
    
    # Apply the boundary conditions
    u[0, :] = u_top
    u[-1, :] = u_bottom
    
    
    u[:, 0] = u_left
    u[:, -1] = u_right
    
    u[Nx//2 -2  : Ny//2 +2,0] = 10

    # Solve the Poisson equation
    for it in range(100):  # Number of iterations
        u_new = np.copy(u)
        for i in range(1, Ny - 1):
            for j in range(1, Nx - 1):
                u_new[i, j] = (u[i-1, j] + u[i+1, j] + u[i, j-1] + u[i, j+1]) / 4
    
        u = np.copy(u_new)
    
    # Calculate the gradients
    grad_x = np.gradient(u, dx, axis=1)
    grad_y = np.gradient(u, dy, axis=0)
    print(grad_x)
    
    # Plot the gradients
    fig, ax = plt.subplots(2, 1, figsize=(6, 8))
    im1 = ax[0].imshow(grad_x, cmap='coolwarm')
    ax[0].set_title('Gradient in X-direction')
    ax[0].set_xlabel('x')
    ax[0].set_ylabel('y')
    plt.colorbar(im1, ax=ax[0])
    
    im2 = ax[1].imshow(grad_y, cmap='coolwarm')
    ax[1].set_title('Gradient in Y-direction')
    ax[1].set_xlabel('x')
    ax[1].set_ylabel('y')
    plt.colorbar(im2, ax=ax[1])
    
    plt.tight_layout()
    plt.show()
    return grad_x,grad_y