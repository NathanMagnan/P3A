### Importations
import orbit_modelling_library as orbit
from mpi4py import MPI

# communicator and local rank
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# parameters
time_isot_start = "2050-01-01T00:00:00"
time_isot_end = "2050-02-01T00:00:00"
n_points_t = 10**4
lambda_start = 10**(11)
lambda_end = 10**(12)
n_points_lambda = 5
alpha_start = 10**(-5)
alpha_end = 10**(-4)
n_points_alpha = 5
target_files = "Test_data_generation"
Bodies = [i + 1 for i in range(10)] 
radiation = False
solar_wind = False
reflectivity_sat = 0.5
radius_sat = 0.63
mass_sat = 100
plasma_speed = 450000

n_proc = n_points_lambda * n_points_alpha + 1

if (rank != n_proc - 1):
    Lambda = [lambda_start + (lambda_end - lambda_start) * i / (n_points_lambda - 1) for i in range(n_points_lambda)]
    Alpha = [alpha_start + (alpha_end - alpha_start) * i / (n_points_alpha - 1) for i in range(n_points_alpha)]
    lambdaa = Lambda[int(rank // n_points_alpha)]
    alpha = Alpha[int(rank % n_points_lambda)]
    
    # computation of the observations
    target_file = target_files + "_lambda_=_" + str(lambdaa) + "_alpha_=_" + str(alpha)
    observation = orbit.Orbit_model(time_isot_start, time_isot_end, n_points_t, target_file, Bodies, "Yukawa", alpha, lambdaa, radiation, solar_wind, reflectivity_sat, radius_sat, mass_sat, plasma_speed)
    observation.compute()
else:
    # computation of the model
    target_file = target_files + "_model"
    model = orbit.Orbit_model(time_isot_start, time_isot_end, n_points_t, target_file, Bodies, "Newton", 0, 1, radiation, solar_wind, reflectivity_sat, radius_sat, mass_sat, plasma_speed)
    model.compute()