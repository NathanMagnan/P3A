## Importations
import math as m
import numpy as np
import scipy.constants as constants
import scipy.interpolate as inter
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import integrate
from jplephem.spk import SPK
from astropy.time import Time, TimeDelta
kernel = SPK.open('C:/Users/Nathan/Documents/Pyzo_Workspace/de430.bsp') # A adapter à son ordinateur

## 
class Satellite:
    def __init__(self, epsilon, rayon, masse):
        self.epsilon = epsilon # Reflectivitée
        self.rayon = rayon
        self.masse = masse

## Définition des forces
class Force():
    def __init__(self):
        self.xpp = " A déterminer"
        self.ypp = " A déterminer"
        self.zpp = " A déterminer"
    
    def calculate(self, V, t):
        return("Error : définir calulate pour cette force")

class Gravitational_Force(Force):
    def __init__(self, n_body):
        self.mu = [2.2032 * 10**(13), 3.24859 * 10**(14), 3.986004418 * 10**(14), 4.282837 * 10**(13), 1.26686534 * 10**(17), 3.7931187 * 10**(16), 5.793939 * 10**(15), 6.836529 * 10**(15), 8.71 * 10**(11), 1.32712440018 * 10**(20)] # m^3 / s^2
        self.n_body = n_body
        
    def calculate(self, V, t):
        t_day = t / 86400 # day
        x, y, z, xp, yp , zp = V
        position_attractor = kernel[0, self.n_body].compute(t_day)
        x_attractor, y_attractor, z_attractor = position_attractor * 1000 # m
        
        d_attractor = m.sqrt((x - x_attractor)**2 + (y - y_attractor)**2 + (z - z_attractor)**2)
        
        self.xpp = - self.mu[self.n_body - 1] * (x - x_attractor) / d_attractor**3
        self.ypp = - self.mu[self.n_body - 1] * (y - y_attractor) / d_attractor**3
        self.zpp = - self.mu[self.n_body - 1] * (z - z_attractor) / d_attractor**3
        return(None)

class Force_Radiation(Force):

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
        F = m.pi * self.sat.rayon**2 * (1 + self.sat.epsilon) * P
        
        if (x - x_sun) > 0 :
            theta = m.atan((y - y_sun) / (x - x_sun))
        elif (x - x_sun) < 0 :
            theta = m.atan((y - y_sun) / (x - x_sun)) + m.pi
        elif (y - y_sun) >= 0 :
            theta = m.pi / 2
        else :
            theta = - m.pi /2
        phi = m.atan((z - z_sun)/(m.sqrt((x - x_sun)**2 + (y - y_sun)**2)))
        
        self.zpp = F * m.sin(phi) / self.sat.masse
        self.ypp = F * m.sin(theta) * m.cos(phi) / self.sat.masse
        self.xpp = F * m.cos(theta) * m.cos(phi) / self.sat.masse
        return (None)

class Force_Vent(Force):

    def __init__(self, vitesse, satellite):
        self.n0 = 1.3 * 10**36   ## particules ejectees par seconde
        self.mp = 1.672 * 10**(-27)   ## masse du proton (ejection de plasma de H+)
        self.v = vitesse ## vitesse du plasma ejecte
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
        F = 2 * m.pi * self.sat.rayon**2 * flux
        
        self.xpp = F * relative_velocity_x / (self.sat.masse * relative_velocity)
        self.ypp = F * relative_velocity_y / (self.sat.masse * relative_velocity)
        self.zpp = F * relative_velocity_z / (self.sat.masse * relative_velocity)
        return (None)

## Propagateur
def propagateur(V, t, Forces, useless_arg):
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

## Tests
## Test de propagation keplerienne
# paramètres de simulation
mu = 1.32712440018 * 10**(20)
a = 1.5 * 10**11
e = 0
p = a * (1 - e**2)
n = p / (1 + e)
b = a * m.sqrt(1 - e**2)
duration = 3 * m.sqrt((2 * m.pi)**2 * a**3 / mu)
k = 6 # ordre de précision

# définition des forces
'''on redéfinit localement une classe de force gravitationelle avec un attracteur fixe''' 
class Gravitational_Force_Test_1(Force):
    def __init__(self):
        self.mu = mu
        return
        
    def calculate(self, V, t):
        t_day = t / 86400 # day
        x, y, z, xp, yp , zp = V
        position_attractor = (0, 0, 0)
        x_attractor, y_attractor, z_attractor = position_attractor
        
        d_attractor = m.sqrt((x - x_attractor)**2 + (y - y_attractor)**2 + (z - z_attractor)**2)
        
        self.xpp = - self.mu * (x - x_attractor) / d_attractor**3
        self.ypp = - self.mu * (y - y_attractor) / d_attractor**3
        self.zpp = - self.mu * (z - z_attractor) / d_attractor**3
        return(None)
force = Gravitational_Force_Test_1()
Forces = [force]

# conditions initiales
x0, y0, z0, xp0, yp0, zp0 = n, 0, 0, 0, m.sqrt(mu / p) * (1 + e), 0
V0 = [x0, y0, z0, xp0, yp0, zp0]

# propagation d'orbite
t = np.linspace(0, duration, num = 10**k)
solution = integrate.odeint(propagateur, V0, t, args = (Forces, 2), printmessg = True, rtol = 10**(-11), hmax = 1, h0 = 1)
X, Y, Z = solution[:, 0], solution[:, 1], solution[:, 2]

# calcul de nu(t)
T = []
Nu = np.linspace(0, 2 * m.pi, 10**k)
for nu in Nu:
    E = 2 * m.atan(m.sqrt(1 - e) * m.tan(nu / 2) / m.sqrt(1 + e))
    M = E - e * m.sin(E)
    if M < 0:
        M += 2 * m.pi
    T.append(M * m.sqrt(a**3 / mu))
nu_de_t = inter.interp1d(T, Nu, bounds_error = False, fill_value = 'extrapolate')

# calcul de l'orbite théorique
X_theorique = []
Y_theorique = []
for i in range(len(t)):
    time = t[i]
    nu = nu_de_t(time)
    r = a * (1 - e**2) / (1 + e * m.cos(nu))
    X_theorique.append(r * m.cos(nu))
    Y_theorique.append(r * m.sin(nu))

# dessin
plt.figure(1)
plt.plot(X, Y, 'r')
plt.plot(X_theorique, Y_theorique, 'g')
plt.show()

## Test de propagation sous pression de radiation
# Conditions initiales
x0, y0, z0, xp0, yp0, zp0 = constants.astronomical_unit, 0, 0, 0, 0, 0
V0 = [x0, y0, z0, xp0, yp0, zp0]
duration = 3 * 86400 * 365
k = 6 # ordre de précision

# Définition des forces
'''on redéfinit localement une classe de force de pression radiative avec une source unique fixe'''
class Force_Radiation_Test(Force):

    def __init__ (self, satellite):
        self.L = 1367 * 4 * m.pi * constants.astronomical_unit**2 # W
        self.sat = satellite

    def calculate(self,V,t):
        x, y, z, xp, yp, zp = V
        r = m.sqrt(x**2 + y**2 + z**2)
        
        flux = self.L / (4 * m.pi * r**2) # W.m^-2
        P = flux / constants.c
        F = m.pi * self.sat.rayon**2 * (1 + self.sat.epsilon) * P
        
        if x > 0 :
            theta = m.atan(y / x)
        elif x < 0 :
            theta = m.atan(y / x) + m.pi
        elif y >= 0 :
            theta = m.pi / 2
        else :
            theta = - m.pi /2
        phi = m.atan(z/(m.sqrt(x**2 + y**2)))
        
        self.zpp = F * m.sin(phi) / self.sat.masse
        self.ypp = F * m.sin(theta) * m.cos(phi) / self.sat.masse
        self.xpp = F * m.cos(theta) * m.cos(phi) / self.sat.masse
        return (None)
satellite = Satellite(0.5, 1, 100)
Forces = [Force_Radiation_Test(satellite)]

# Propagation d'orbite
t = np.linspace(0, duration, num = 10**k)
solution = integrate.odeint(propagateur, V0, t, args = (Forces, 2), printmessg = True, rtol = 10**(-11), hmax = 1, h0 = 1)
X, Y, Z = solution[:, 0], solution[:, 1], solution[:, 2]

# Dessin
plt.figure(2)
plt.plot(X, Y, 'r')
plt.show()

## Test de propagation sous vent solaire
# Conditions initiales
x0, y0, z0, xp0, yp0, zp0 = constants.astronomical_unit, 0, 0, 0, 0, 0
V0 = [x0, y0, z0, xp0, yp0, zp0]
duration = 365 * 86400
k = 6 # ordre de précision

# Définition des forces
'''on redéfinit localement une classe de force du vent solaire avec une source unique fixe'''
class Force_Vent_Test(Force):

    def __init__(self, vitesse, satellite):
        self.n0 = 1.3 * 10**36   ## particules ejectees par seconde
        self.mp = 1.672 * 10**(-27)   ## masse du proton (ejection de plasma de H+)
        self.v = vitesse ## vitesse du plasma ejecte
        self.sat = satellite

    def calculate(self,V,t):
        x, y, z, xp, yp, zp = V
        r = m.sqrt(x**2 + y**2 + z**2)
        flux = self.mp * self.v * self.n0 / (4 * m.pi * r**2)
        F = 2 * m.pi * self.sat.rayon**2 * flux
        
        if x > 0 :
            theta = m.atan(y/x)
        elif x < 0 :
            theta = m.atan(y / x) + m.pi
        elif y >= 0 :
            theta = m.pi / 2
        else :
            theta = - m.pi /2
        phi = m.atan(z/(m.sqrt(x**2 + y**2)))
        
        self.zpp = F * m.sin(phi) / self.sat.masse
        self.ypp = F * m.sin(theta) * m.cos(phi) / self.sat.masse
        self.xpp = F * m.cos(theta) * m.cos(phi) / self.sat.masse
        return (None)
satellite = Satellite(0.5, 1, 100)
vitesse = 450000
Forces = [Force_Vent_Test(vitesse, satellite)]

# Propagation d'orbite
t = np.linspace(0, duration, num = 10**k)
solution = integrate.odeint(propagateur, V0, t, args = (Forces, 2), printmessg = True, rtol = 10**(-11), hmax = 1, h0 = 1)
X, Y, Z = solution[:, 0], solution[:, 1], solution[:, 2]

# Dessin
plt.figure(3)
plt.plot(X, Y, 'b')
plt.show()

## Test de propagation sous le potentiel de Yukawa

## Test des niveaux des forces
# Définition des forces
'''on redéfinit localement une classe de force gravitationelle avec un attracteur au niveau de la Terre''' 
class Gravitational_Force_Test_2(Force):
    def __init__(self, n_body):
        self.mu = [2.2032 * 10**(13), 3.24859 * 10**(14), 3.986004418 * 10**(14), 4.282837 * 10**(13), 1.26686534 * 10**(17), 3.7931187 * 10**(16), 5.793939 * 10**(15), 6.836529 * 10**(15), 8.71 * 10**(11), 1.32712440018 * 10**(20)] # m^3 / s^2
        self.n_body = n_body
        
    def calculate(self, V, t):
        t_day = t / 86400 # day
        x, y, z, xp, yp , zp = V
        position_attractor = kernel[0, self.n_body].compute(t_day)
        position_attractor -= kernel[0, 3].compute(t_day)
        #position_attractor -= kernel[3, 399].compute(t_day)
        x_attractor, y_attractor, z_attractor = position_attractor * 1000 # m
        
        d_attractor = m.sqrt((x - x_attractor)**2 + (y - y_attractor)**2 + (z - z_attractor)**2)
        
        self.xpp = - self.mu[self.n_body - 1] * (x - x_attractor) / d_attractor**3
        self.ypp = - self.mu[self.n_body - 1] * (y - y_attractor) / d_attractor**3
        self.zpp = - self.mu[self.n_body - 1] * (z - z_attractor) / d_attractor**3
        return(None)
force_earth = Gravitational_Force_Test_2(3)
force_sun = Gravitational_Force_Test_2(10)
Forces = [force_earth, force_sun]

# Définition des temps
time_isot = "2050-01-01T00:00:00"
Times_isot = [time_isot]
T = Time(Times_isot, format = 'isot', scale = 'utc')
Times_jd = T.jd
t = Times_jd[0] * 86400

# Position de référence dans le référentiel de la Terre (399)
x0, y0, z0, xp0, yp0, zp0 = 0, 0, 0, 0, 0, 0
V0 = [x0, y0, z0, xp0, yp0, zp0]

# Calcul des forces
H = []
Earth = []
Sun = []
for h in range(6400, 100000, 100):
    x, y, z, vx, vy, vz = V0
    x = 1000 * h
    V = [x, y, z, vx, vy, vz]
    
    Forces[0].calculate(V, t)
    xpp, ypp, zpp = Forces[0].xpp, Forces[0].ypp, Forces[0].zpp 
    a_earth = m.sqrt(xpp**2 + ypp**2 + zpp**2)
    Forces[1].calculate(V, t)
    xpp, ypp, zpp = Forces[1].xpp, Forces[1].ypp, Forces[1].zpp
    a_sun = m.sqrt(xpp**2 + ypp**2 + zpp**2)
    
    H.append(h)
    Earth.append(a_earth)
    Sun.append(a_sun)

# Dessin
plt.figure(4)
plt.xscale("log")
plt.yscale("log")
plt.plot(H, Earth, 'b')
plt.plot(H, Sun, 'r')
plt.show()

## Calculs
## Influence de Jupyter
# définition des temps
time_isot_start = "2050-01-01T00:00:00"
time_isot_end = "2050-01-21T00:00:00"
Times_isot = [time_isot_start, time_isot_end]
T = Time(Times_isot, format = 'isot', scale = 'utc')
Times_jd = T.jd
time_jd_start, time_jd_end = Times_jd
time_jd_start *= 86400
time_jd_end *= 86400
n_points = 10**6

# conditions initiales
position0, velocity0 = kernel[0, 9].compute_and_differentiate(time_jd_start / 86400)
x0, y0, z0 = position0 * 1000 # m
x0 += 10**(11) # on s'éloigne de Pluton d'une unité astronomique
xp0, yp0, zp0 = velocity0 * 1000 / 86400 # m / s
V0 = [x0, y0, z0, xp0, yp0, zp0]

# définition des forces avec tous les corps
Forces_sun_and_jupyter = []
for n_body in range(1,11):
    if ((n_body == 10) or (n_body == 5)):
        force = Gravitational_Force(n_body)
        Forces_sun_and_jupyter.append(force)
# définition des forces avec juste le soleil
Forces_sun = []
for n_body in range(1,11):
    if (n_body == 10):
        force = Gravitational_Force(n_body)
        Forces_sun.append(force)

# propagation d'orbite
t = np.linspace(time_jd_start, time_jd_end, n_points)
print("Start")

V01 = V0[:]
solution = integrate.odeint(propagateur, V01, t, args = (Forces_sun_and_jupyter, 2), printmessg = True, rtol = 10**(-11), hmax = 1, h0 = 1)
X_sun_and_jupyter, Y_sun_and_jupyter, Z_sun_and_jupyter = solution[:, 0], solution[:, 1], solution[:, 2]
print("OK sun_and_jupyter")

V02 = V0[:]
solution_2 = integrate.odeint(propagateur, V02, t, args = (Forces_sun, 2), printmessg = True, rtol = 10**(-11), hmax = 1, h0 = 1)
X_sun, Y_sun, Z_sun = solution_2[:, 0], solution_2[:, 1], solution_2[:, 2]
print("OK sans")

# dessin
plt.figure(5)
plt.plot(X_sun_and_jupyter, Y_sun_and_jupyter, 'g')
plt.plot(X_sun, Y_sun, 'r')
plt.show()
print(m.sqrt((X_sun_and_jupyter[-1] - X_sun[-1])**2 + (Y_sun_and_jupyter[-1] - Y_sun[-1])**2 + (Z_sun_and_jupyter[-1] - Z_sun[-1])**2))

## Influence des autres planètes
# définition des temps
time_isot_start = "2050-01-01T00:00:00"
time_isot_end = "2053-01-01T00:00:00"
Times_isot = [time_isot_start, time_isot_end]
T = Time(Times_isot, format = 'isot', scale = 'utc')
Times_jd = T.jd
time_jd_start, time_jd_end = Times_jd
time_jd_start *= 86400
time_jd_end *= 86400
n_points = 10**6

# conditions initiales
position0, velocity0 = kernel[0, 9].compute_and_differentiate(time_jd_start / 86400)
x0, y0, z0 = position0 * 1000 # m
x0 += 10**(11) # on s'éloigne de Pluton d'une unité astronomique
xp0, yp0, zp0 = velocity0 * 1000 / 86400 # m / s
V0 = [x0, y0, z0, xp0, yp0, zp0]

# définition des forces avec toutes les planètes
Forces_all_planets = []
for n_body in range(1,11):
    if (n_body != 3):
        force = Gravitational_Force(n_body)
        Forces_all_planets.append(force)

# propagation d'orbite
t = np.linspace(time_jd_start, time_jd_end, n_points)
print("Start")

V03 = V0[:]
solution = integrate.odeint(propagateur, V03, t, args = (Forces_all_planets, 2), printmessg = True, rtol = 10**(-11), hmax = 1, h0 = 1)
X_all_planets, Y_all_planets, Z_all_planets = solution[:, 0], solution[:, 1], solution[:, 2]
print("OK all_planets")

# dessin
plt.figure(6)
plt.plot(X_all_planets, Y_all_planets, 'g')
plt.plot(X_sun_and_jupyter, Y_sun_and_jupyter, 'r')
plt.show()
print(m.sqrt((X_all_planets[-1] - X_sun_and_jupyter[-1])**2 + (Y_all_planets[-1] - Y_sun_and_jupyter[-1])**2 + (Z_all_planets[-1] - Z_sun_and_jupyter[-1])**2))

## Influence des autres forces
# définition des temps
time_isot_start = "2050-01-01T00:00:00"
time_isot_end = "2053-01-01T00:00:00"
Times_isot = [time_isot_start, time_isot_end]
T = Time(Times_isot, format = 'isot', scale = 'utc')
Times_jd = T.jd
time_jd_start, time_jd_end = Times_jd
time_jd_start *= 86400
time_jd_end *= 86400
n_points = 10**6

# définition du satellite
satellite = Satellite(0.5, 0.63, 1000)

# conditions initiales
position0, velocity0 = kernel[0, 9].compute_and_differentiate(time_jd_start / 86400)
x0, y0, z0 = position0 * 1000 # m
x0 += 10**(11) # on s'éloigne de Pluton d'une unité astronomique
xp0, yp0, zp0 = velocity0 * 1000 / 86400 # m / s
V0 = [x0, y0, z0, xp0, yp0, zp0]

# définition des forces
Forces_all_forces = []
for n_body in range(1,11):
    if (n_body != 3):
        force = Gravitational_Force(n_body)
        Forces_all_forces.append(force)
force = Force_Radiation(satellite)
Forces_all_forces.append(force)
force = Force_Vent(450000, satellite)
Forces_all_forces.append(force)

# propagation d'orbite
t = np.linspace(time_jd_start, time_jd_end, n_points)
print("Start")

V04 = V0[:]
solution = integrate.odeint(propagateur, V03, t, args = (Forces_all_forces, 2), printmessg = True, rtol = 10**(-11), hmax = 1, h0 = 1)
X_all_forces, Y_all_forces, Z_all_forces = solution[:, 0], solution[:, 1], solution[:, 2]
print("OK all_forces")

# dessin
plt.figure(7)
plt.plot(X_all_forces, Y_all_forces, 'g')
plt.plot(X_all_planets, Y_all_planets, 'r')
plt.show()
print(m.sqrt((X_all_forces[-1] - X_all_planets[-1])**2 + (Y_all_forces[-1] - Y_all_planets[-1])**2 + (Z_all_forces[-1] - Z_all_planets[-1])**2))

## Différence entre le pootentiel de Yukawa et le potentiel de Newton
