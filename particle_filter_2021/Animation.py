import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.animation as animation
import copy

from Simulator import Simulator
from SimulationConst import SimulationConst

"""
This script is ment for making an animation of the particle filter simulation

Do not edit this
"""
class Animation(animation.TimedAnimation):
    def __init__(self, N_particles, x_r, y_r, phi, kappa, state, sense, simConst: SimulationConst, fps=120):
        self.interval = 1000/fps

        fig = plt.figure()
        ax = fig.add_subplot(111)

        self.N_particles = N_particles
        self.x_r   = x_r
        self.y_r   = y_r
        self.phi   = phi
        self.kappa = kappa

        self.K     = state.shape[0]
        self.state = state

        self.sense = sense

        self.simConst_contour = simConst.contour
        self.contour = copy.deepcopy(simConst.contour)

        self.contour[7,0] += self.state[0,3]
        self.contour[8,0] += self.state[0,3]

        self.robot_points = 0.05*np.array([[-1,-1],[-1,1],[0.5,1],[1,0],[0.5,-1],[-1,-1]])

        self.contour_plot     = Line2D([], [], color='black', linestyle='--')
        self.robot_area       = Line2D([], [], color='green')
        self.distance_line    = Line2D([], [], color='red')
        self.particle_circles = Line2D([], [], color='red', marker='o', linestyle='')
        self.particle_walls   = Line2D([], [], color='red', alpha=0.75)
        
        ax.add_line(self.contour_plot)
        ax.add_line(self.robot_area)
        ax.add_line(self.distance_line)
        ax.add_line(self.particle_circles)
        ax.add_line(self.particle_walls)

        spacing = 0.6
        ax.set_xlim(min(self.contour[:,0])-spacing, max(self.contour[:,0])+spacing)
        ax.set_ylim(min(self.contour[:,1])-spacing, max(self.contour[:,1])+spacing)
        ax.set_aspect('equal')

        animation.TimedAnimation.__init__(self, fig, interval=self.interval, blit=False)

    def _draw_frame(self, framedata):
        k = framedata

        self.robot_area.set_data(np.cos(self.state[k,2])*self.robot_points[:,0]-np.sin(self.state[k,2])*self.robot_points[:,1]+self.state[k,0],
                                 np.sin(self.state[k,2])*self.robot_points[:,0]+np.cos(self.state[k,2])*self.robot_points[:,1]+self.state[k,1])

        self.particle_circles.set_data(self.x_r[k,:], self.y_r[k,:])

        sort_kappa = np.sort(self.kappa[k,:]) + self.simConst_contour[7,0]
        self.particle_walls.set(xdata=np.repeat(sort_kappa, 2))

        if k > 1:
            self.distance_line.set_data(np.cos(self.state[k,2])*np.hstack((0.0, self.sense[k,:]))+self.state[k,0],
                                        np.sin(self.state[k,2])*np.hstack((0.0, self.sense[k,:]))+self.state[k,1])

        self._drawn_artists = [self.robot_area, self.particle_circles, self.particle_walls, self.distance_line]

    def new_frame_seq(self):
        return iter(range(self.K))

    def _init_draw(self):
        lines = [self.robot_area, self.particle_circles, self.distance_line]
        for l in lines:
            l.set_data([], [])

        #preset ydata of particles walls
        self.particle_walls.set(ydata=np.tile([self.contour[7,1], self.contour[8,1]], self.N_particles))

        #contour plot
        xdata = np.hstack((self.contour[:,0], self.contour[0,0]))
        ydata = np.hstack((self.contour[:,1], self.contour[0,1]))
        self.contour_plot.set_data(xdata, ydata)
