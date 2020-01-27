## Importations
import orbit_modelling_library as orbit
import math as m
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D
rc('text', usetex = True)

## Fonction chi2
def chi_2(observation, model, delta_radar, delta_accelerometer, delta_initial_position, delta_initial_speed):
    if len(observation.t) != len(model.t):
        print("Error : the observation and the model are not sampled at the same times")
        return
    for i in range(len(observation.t)):
        if observation.t[i] != model.t[i]:
            print("Error : the observation and the model are not sampled at the same times")
            return
    
    chi_2 = 0
    delta_position = delta_initial_position
    delta_speed = delta_initial_speed
    T = observation.t
    t0 = T[0]
    time_step = T[1] - T[0]
    
    for i in range(len(T)):
        X_obs, X_mod = observation.X[i], model.X[i]
        Y_obs, Y_mod = observation.Y[i], model.Y[i]
        Z_obs, Z_mod = observation.Z[i], model.Z[i]
        
        delta_position += time_step * delta_speed
        delta_speed += time_step * delta_accelerometer
        
        chi_2 += ((X_obs - X_mod)**2 + (Y_obs - Y_mod)**2 + (Z_obs - Z_mod)**2) / (delta_radar**2 + delta_position**2)
    return chi_2

## Data generation
def observations_data_generation(time_isot_start, time_isot_end, n_points_t, 
                                      lambda_start, lambda_end, n_points_lambda, 
                                      alpha_start, alpha_end, n_points_alpha, 
                                      target_files, 
                                      Bodies = [i + 1 for i in range(10)], radiation = False, solar_wind = False, 
                                      reflectivity_sat = 0.5, radius_sat = 0.63, mass_sat = 100, plasma_speed = 450000):
    Lambda = [lambda_start + (lambda_end - lambda_start) * i / (n_points_lambda - 1) for i in range(n_points_lambda)]
    Alpha = [alpha_start + (alpha_end - alpha_start) * i / (n_points_alpha - 1) for i in range(n_points_alpha)]
    
    N = n_points_lambda*n_points_alpha
    print(N)
    
    i = 0
    print(100 * i / (N - 1))
    for lambdaa in Lambda:
        for alpha in Alpha:
            i += 1
            target_file = target_files + "_lambda_=_" + str(lambdaa) + "_alpha_=_" + str(alpha)
            observation = orbit.Orbit_model(time_isot_start, time_isot_end, n_points_t, target_file, Bodies, "Yukawa", alpha, lambdaa, radiation, solar_wind, reflectivity_sat, radius_sat, mass_sat, plasma_speed)
            observation.compute()
            print(100 * i / N) # gives the pourcentage of the computation that has already been done

def model_data_generation(time_isot_start, time_isot_end, n_points_t, target_files, 
                          Bodies = [i + 1 for i in range(10)], type = "Yukawa", alpha = 10**(-4), lambdaa = 10**12, 
                          radiation = False, solar_wind = False, 
                          reflectivity_sat = 0.5, radius_sat = 0.63, mass_sat = 100, plasma_speed = 450000):
    
    target_file = target_files + "_model"
    model = orbit.Orbit_model(time_isot_start, time_isot_end, n_points_t, target_file, Bodies, type, alpha, lambdaa, radiation, solar_wind, reflectivity_sat, radius_sat, mass_sat, plasma_speed)
    model.compute()

def observations_data_generation(time_isot_start, time_isot_end, n_points_t, 
                                      lambda_start, lambda_end, n_points_lambda, 
                                      alpha_start, alpha_end, n_points_alpha, 
                                      target_files, 
                                      Bodies = [i + 1 for i in range(10)], radiation = False, solar_wind = False, 
                                      reflectivity_sat = 0.5, radius_sat = 0.63, mass_sat = 100, plasma_speed = 450000):
    Lambda = [lambda_start + (lambda_end - lambda_start) * i / (n_points_lambda - 1) for i in range(n_points_lambda)]
    Alpha = [alpha_start + (alpha_end - alpha_start) * i / (n_points_alpha - 1) for i in range(n_points_alpha)]
    
    N = n_points_lambda*n_points_alpha
    print(N)
    
    i = 0
    print(100 * i / (N - 1))
    for lambdaa in Lambda:
        for alpha in Alpha:
            i += 1
            target_file = target_files + "_lambda_=_" + str(lambdaa) + "_alpha_=_" + str(alpha)
            observation = orbit.Orbit_model(time_isot_start, time_isot_end, n_points_t, target_file, Bodies, "Yukawa", alpha, lambdaa, radiation, solar_wind, reflectivity_sat, radius_sat, mass_sat, plasma_speed)
            observation.compute()
            print(100 * i / N) # gives the pourcentage of the computation that has already been done

observations_data_generation("2050-01-01T00:00:00", "2050-02-01T00:00:00", 10**4, 
                            10**(-12), 10**(-11), 5, 
                            10**(-5), 10**(-4), 5, 
                            "Test_data_generation")