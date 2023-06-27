import numpy as np

class PostParticles:
    def __init__(self, N_particles=1000):
        self.N_particles = N_particles

        self.x_r = np.empty(self.N_particles)
        self.y_r = np.empty(self.N_particles)

        self.phi = np.empty(self.N_particles)
        self.kappa = np.empty(self.N_particles)