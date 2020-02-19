## Importations
import math as m
import numpy as np
import scipy.constants as constants
import scipy.interpolate as inter
import matplotlib.pyplot as plt
from scipy import integrate
from jplephem.spk import SPK
from astropy.time import Time, TimeDelta
kernel = SPK.open(r'C:/Users/Nathan/Documents/W - Pyzo_Workspace/de430.bsp') # To be adapted to one's computer

## Code
class Satellite:
    
    def __init__(self, reflectivity, radius, mass):
        self.reflectivity = reflectivity # epsilon
        self.radius = radius
        self.mass = mass
 
        
class Force():
    
    def __init__(self):
        self.xpp = " To be determined"
        self.ypp = " TO be determined"
        self.zpp = " To be determined"
    
    def calculate(self, V, t):
        return("Error : define calulate for this force")


class Gravitational_Force(Force):
    
    def __init__(self, n_body):
        self.Mu = [2.2032 * 10**(13), 3.24859 * 10**(14), 3.986004418 * 10**(14), 4.282837 * 10**(13), 1.26686534 * 10**(17), 3.7931187 * 10**(16), 5.793939 * 10**(15), 6.836529 * 10**(15), 8.71 * 10**(11), 1.32712440018 * 10**(20)] # m^3 / s^2
        self.n_body = n_body
        self.mu = self.Mu[n_body - 1]
        
    def calculate(self, V, t):
        t_day = t / 86400 # day
        x, y, z, xp, yp , zp = V
        position_attractor = kernel[0, self.n_body].compute(t_day)
        x_attractor, y_attractor, z_attractor = position_attractor * 1000 # m
        
        d_attractor = m.sqrt((x - x_attractor)**2 + (y - y_attractor)**2 + (z - z_attractor)**2)
        
        self.xpp = - self.mu * (x - x_attractor) / d_attractor**3
        self.ypp = - self.mu * (y - y_attractor) / d_attractor**3
        self.zpp = - self.mu * (z - z_attractor) / d_attractor**3
        return(None)


class Yukawa_Force(Force):
    
    def __init__(self, n_body, alpha, lambdaa):
        self.Mu = [2.2032 * 10**(13), 3.24859 * 10**(14), 3.986004418 * 10**(14), 4.282837 * 10**(13), 1.26686534 * 10**(17), 3.7931187 * 10**(16), 5.793939 * 10**(15), 6.836529 * 10**(15), 8.71 * 10**(11), 1.32712440018 * 10**(20)] # m^3 / s^2
        self.n_body = n_body
        self.mu = self.Mu[n_body - 1]
        
        self.alpha = alpha 
        self.lambdaa = lambdaa
    
    def calculate(self, V, t):
        t_day = t / 86400 # day
        
        x, y, z, xp, yp , zp = V
        position_attractor = kernel[0, self.n_body].compute(t_day)
        x_attractor, y_attractor, z_attractor = position_attractor * 1000 # m
        
        d_attractor = m.sqrt((x - x_attractor)**2 + (y - y_attractor)**2 + (z - z_attractor)**2)
        F = - self.mu * d_attractor**(-2) * (1 + self.alpha * m.exp(-d_attractor / self.lambdaa) * (1 + d_attractor / self.lambdaa))
        
        theta = m.atan2((y - y_attractor), (x - x_attractor))
        phi = m.atan((z - z_attractor) / (m.sqrt((x - x_attractor)**2 + (y - y_attractor)**2)))
        
        self.xpp = F * m.cos(theta) * m.cos(phi)
        self.ypp = F * m.sin(theta) * m.cos(phi) 
        self.zpp = F * m.sin(phi)
        return(None)


class Radiation_Force(Force):

    def __init__ (self, satellite):
        self.L = 1367 * 4 * m.pi * constants.astronomical_unit**2 # W
        self.sat = satellite

    def calculate(self,V,t):
        t_day = t / 86400 # day
        
        x, y, z, xp, yp, zp = V
        
        position_sun = kernel[0, 9].compute(t_day)
        x_sun, y_sun, z_sun = position_sun * 1000 # m
        
        d_sun = m.sqrt((x - x_sun)**2 + (y - y_sun)**2 + (z - z_sun)**2)
        
        flux = self.L / (4 * m.pi * d_sun**2) # W.m^-2
        P = flux / constants.c
        F = m.pi * self.sat.radius**2 * (1 + self.sat.reflectivity) * P
        
        theta = m.atan2(y - y_sun, x - x_sun)
        phi = m.atan((z - z_sun)/(m.sqrt((x - x_sun)**2 + (y - y_sun)**2)))
        
        self.zpp = F * m.sin(phi) / self.sat.mass
        self.ypp = F * m.sin(theta) * m.cos(phi) / self.sat.mass
        self.xpp = F * m.cos(theta) * m.cos(phi) / self.sat.mass
        return (None)


class Solar_Wind_Force(Force):

    def __init__(self, plasma_speed, satellite):
        self.n0 = 1.3 * 10**36   ## particles ejected per seconde
        self.mp = 1.672 * 10**(-27)   ## mass of the proton (with assume a plasma of H+)
        self.v = plasma_speed ## speed of the ejected plasma
        self.sat = satellite

    def calculate(self,V,t):
        t_day = t / 86400 # day
        
        x, y, z, xp, yp, zp = V
        
        position_sun = kernel[0, 9].compute(t_day)
        x_sun, y_sun, z_sun = position_sun * 1000 # m
        d_sun = m.sqrt((x - x_sun)**2 + (y - y_sun)**2 + (z - z_sun)**2)
        
        relative_velocity_x = self.v * (x - x_sun) / d_sun - xp
        relative_velocity_y = self.v * (y - y_sun) / d_sun - yp
        relative_velocity_z = self.v * (z - z_sun) / d_sun - zp
        relative_velocity = m.sqrt(relative_velocity_x**2 + relative_velocity_y**2 + relative_velocity_z**2)
        
        flux = self.mp * relative_velocity * self.n0 / (4 * m.pi * d_sun**2)
        F = 2 * m.pi * self.sat.radius**2 * flux
        
        self.xpp = F * relative_velocity_x / (self.sat.mass * relative_velocity)
        self.ypp = F * relative_velocity_y / (self.sat.mass * relative_velocity)
        self.zpp = F * relative_velocity_z / (self.sat.mass * relative_velocity)
        return (None)


def propagator(V, t, Forces, useless_arg):
    x, y, z, xp, yp, zp = V
    xpp = 0
    ypp = 0
    zpp = 0
    for f in Forces:
        f.calculate(V, t)
        xpp += f.xpp
        ypp += f.ypp
        zpp += f.zpp
    return [xp, yp, zp, xpp, ypp, zpp]

## Encapsulation

class Orbit_model():
    def __init__(self, time_isot_start, time_isot_end, n_points, target_file, Bodies = [i + 1 for i in range(10)], type = "Yukawa", alpha = 10**(-4), lambdaa = 10**12, radiation = False, solar_wind = False, reflectivity_sat = 0.5, radius_sat = 0.63, mass_sat = 100, plasma_speed = 450000):
        self.sat = Satellite(reflectivity_sat, radius_sat, mass_sat)
        
        self.Forces = []
        if radiation:
            self.Forces.append(Radiation_Force(self.sat))
        if solar_wind:
            self.Forces.append(Solar_Wind_Force(plasma_speed, self.sat))
        for n_body in Bodies:
            if type == "Yukawa":
                self.Forces.append(Yukawa_Force(n_body, alpha, lambdaa))
            else :
                self.Forces.append(Gravitational_Force(n_body))
        
        Times_isot = [time_isot_start, time_isot_end]
        T = Time(Times_isot, format = 'isot', scale = 'utc')
        Times_jd = T.jd
        self.time_jd_start, self.time_jd_end = Times_jd
        
        self.time_jd_start *= 86400
        self.time_jd_end *= 86400
        self.n_points = n_points
        
        self.t = np.linspace(self.time_jd_start, self.time_jd_end, self.n_points)
        self.delta_t = self.t[1] - self.t[0]
        
        position0, velocity0 = kernel[0, 9].compute_and_differentiate(self.time_jd_start / 86400)
        x0, y0, z0 = position0 * 1000 # m
        x0 += 10**(11) # we get away from pluto by one astronomical unit
        xp0, yp0, zp0 = velocity0
        wrong_celerity = m.sqrt(xp0**2 + yp0**2 + zp0**2)
        right_celerity = constants.c * 0.0005
        xp0 = xp0 * right_celerity / wrong_celerity # m /s
        yp0 = yp0 * right_celerity / wrong_celerity # m /s
        zp0 = zp0 * right_celerity / wrong_celerity # m /s
        self.V0 = [x0, y0, z0, xp0, yp0, zp0]
        
        self.X = np.asarray([])
        self.Y = np.asarray([])
        self.Z = np.asarray([])
        
        self.target = target_file
    
    def __str__(self): # To be done
        return "To Be Done : __str__ for class Orbit_model"
    
    def save(self):
        np.savetxt(str(self.target) + "_X", self.X)
        np.savetxt(str(self.target) + "_Y", self.Y)
        np.savetxt(str(self.target) + "_Z", self.Z)
        return
    
    def load(self):
        self.X = np.loadtxt(str(self.target) + "_X")
        self.Y = np.loadtxt(str(self.target) + "_Y")
        self.Z = np.loadtxt(str(self.target) + "_Z")
        if len(self.X) != len(self.t) :
            print("Error : the target files does not match the model start and end date")
            return
        else:
            return
    
    def plot(self, color): # To be tested
        self.load()
        ax = plt.plot(self.X, self.Y, color)
        return ax
    
    def compute(self):
        solution = integrate.odeint(propagator, self.V0, self.t, args = (self.Forces, 2), printmessg = True, rtol = 10**(-11), hmax = 1, h0 = 1)
        self.X, self.Y, self.Z = solution[:, 0], solution[:, 1], solution[:, 2]
        self.save()
        return    