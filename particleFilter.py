import numpy as np
from filterpy.monte_carlo import systematic_resample
from numpy.linalg import norm
from numpy.random import randn
import matplotlib.pyplot as plt
import json
import navpy
from pyproj import Transformer
from scipy.stats import multivariate_normal
from scipy.linalg import inv

def create_particles(geoSpatialTransforms, N):
    particles = np.empty((N, 3))
    lat, lon, alt = geoSpatialTransforms[0][:3]
    x, y, z = latlonalt_to_ecef(lat, lon, alt)
    particles[:, 0] = x
    particles[:, 1] = y
    particles[:, 2] = z
    print(particles)
    return particles

def latlonalt_to_ecef(lat, lon, alt):
    transformer = Transformer.from_crs("epsg:4326", "epsg:4978")
    x, y, z = transformer.transform(lat, lon, alt)
    return x, y, z

def update_ecef_uncertainties(lat, lon, alt, lat_uncertainty, lon_uncertainty, alt_uncertainty):
    a = 6378137.0  # Earth's semi-major axis (m)
    f = 1 / 298.257223563  # Earth's flattening
    e2 = 2 * f - f ** 2  # Earth's first eccentricity squared
    #Not very good
    N_phi = a / np.sqrt(1 - e2 * np.sin(np.radians(lat)) ** 2)
    x_uncertainty = np.abs((N_phi + alt) * np.cos(np.radians(lat)) * np.radians(lon_uncertainty))
    y_uncertainty = np.abs((N_phi * (1 - e2) + alt) * np.sin(np.radians(lat)) * np.radians(lat_uncertainty))
    z_uncertainty = alt_uncertainty

    return x_uncertainty, y_uncertainty, z_uncertainty
def predict(particles, prev_lat, prev_lon, prev_alt, current_lat, current_lon, current_alt, prev_pose, pose, noise_std=1):
    prev_x, prev_y, prev_z = latlonalt_to_ecef(prev_lat, prev_lon, prev_alt)
    current_x, current_y, current_z = latlonalt_to_ecef(current_lat, current_lon, current_alt)
    dist_ecef = np.array([current_x - prev_x, current_y - prev_y, current_z - prev_z])
    noise = np.random.normal(0, noise_std, size=(particles.shape[0], 3))
    particles += dist_ecef  + noise
# def predict(particles, delta_time, prev_lat, prev_lon, prev_alt, heading, avg_walking_speed=1.4, noise_std=20):
#     R = 6371000  # radius of Earth in meters
#     lat_rad = np.radians(prev_lat)
#     lon_rad = np.radians(prev_lon)
    
#     dlat = (avg_walking_speed * delta_time * np.cos(np.radians(heading))) / R
#     dlon = (avg_walking_speed * delta_time * np.sin(np.radians(heading))) / (R * np.cos(lat_rad))
    
#     new_lat = prev_lat + np.degrees(dlat)
#     new_lon = prev_lon + np.degrees(dlon)

#     prev_x, prev_y, prev_z = latlonalt_to_ecef(prev_lat, prev_lon, prev_alt)
#     current_x, current_y, current_z = latlonalt_to_ecef(new_lat, new_lon, prev_alt)
    
#     dist_ecef = np.array([current_x - prev_x, current_y - prev_y, current_z - prev_z])

#     noise = np.random.normal(0, noise_std, size=(particles.shape[0], 3))
#     particles += dist_ecef + noise


def update(particles, weights, measurement_ecef, uncertainties, weight_effect=4000000):
    weights.fill(1.0)

    for i, p in enumerate(particles):
        cov = np.diag(np.array(uncertainties) ** 2)
        weight = multivariate_normal.pdf(p, mean=measurement_ecef, cov=cov)
        weights[i] = weight ** weight_effect
    weights += 1.e-300
    weights /= np.sum(weights)
    return weights

def estimate(particles, weights):
    return np.average(particles, axis=0, weights=weights)

def neff(weights):
    return 1. / np.sum(np.square(weights))

def resample_from_index(particles, weights, indexes):
    particles[:] = particles[indexes]
    weights.resize(len(particles))
    weights.fill(1.0 / len(weights))

def particle_filter(geoSpatialTransforms, N=1000):
    particles = create_particles(geoSpatialTransforms, N)
    weights = np.ones(N) / N
    predictions = []

    for t in range(1, len(geoSpatialTransforms)):
        prev_lat, prev_lon, prev_alt = geoSpatialTransforms[t - 1][:3]
        prev_pose = geoSpatialTransforms[t-1][-1]
        pose = geoSpatialTransforms[t][-1]
        current_lat, current_lon, current_alt = geoSpatialTransforms[t][:3]
        heading = geoSpatialTransforms[t][3]
        #predict(particles, delta_time, prev_lat, prev_lon, prev_alt, heading)
        predict(particles, prev_lat, prev_lon, prev_alt, current_lat, current_lon, current_alt, prev_pose, pose)
        print(f"Step: {t}, Heading: {heading}")
        
        measurement_ecef = latlonalt_to_ecef(current_lat, current_lon, current_alt)

        uncertainties = update_ecef_uncertainties(current_lat, current_lon, current_alt,
                                          geoSpatialTransforms[t][4],
                                          geoSpatialTransforms[t][6],
                                          geoSpatialTransforms[t][3])
        weights = update(particles, weights, measurement_ecef, uncertainties)
        

        print(f"Weights at step {t}: {weights}")
        
        if neff(weights) < N / 2:
            return
            indexes = systematic_resample(weights)
            resample_from_index(particles, weights, indexes)
        
        print(f"Particles at step {t}: {particles}")

        estimate_position = estimate(particles, weights)
        
        print(f"Estimated position at step {t}: {estimate_position}")
        
        predictions.append(estimate_position)

    return predictions

f = open('bad_metadata.json')
metadata = json.load(f)
f2 = open('bad_pathdata.json')
pathdata = json.load(f2)

alldata = pathdata
for key in metadata:
    alldata[key] = metadata[key]
timestamps = alldata["geoSpatialTransformTimes"]
coords = []
for d in alldata["garAnchorCameraWorldTransformsAndGeoSpatialData"]:
    lat = d["geospatialTransform"]['latitude']
    lon = d["geospatialTransform"]['longitude']
    alt = d["geospatialTransform"]['altitude']
    heading = d["geospatialTransform"]['heading']
    lat_uncertainty = d["geospatialTransform"]['positionAccuracy']
    lon_uncertainty = d["geospatialTransform"]['positionAccuracy']
    alt_uncertainty = d["geospatialTransform"]['altitudeAccuracy']
    pose = np.array(d["cameraWorldTransform"]).reshape(4, 4).T
    coords.append([lat, lon, alt, heading, lat_uncertainty, lon_uncertainty, alt_uncertainty, pose])

transform_poses = [np.array(x["cameraWorldTransform"]).reshape(4, 4).T for x in alldata["garAnchorCameraWorldTransformsAndGeoSpatialData"]]
inv(transform_poses[1]) @ transform_poses[0]

predictions = particle_filter(coords)

# ecef_coords = [latlonalt_to_ecef(coord[0], coord[1], coord[2]) for coord in coords]
# x_coords = [coord[0] for coord in ecef_coords]
# y_coords = [coord[1] for coord in ecef_coords]

# plt.plot(x_coords, y_coords, 'go', label='Ground Truth')
# plt.plot([pred[0] for pred in predictions], [pred[1] for pred in predictions], 'bx', label='Predicted')
# plt.legend(loc='lower right')
# plt.xlabel('ECEF X')
# plt.ylabel('ECEF Y')
# plt.title('Particle Filter Predictions in ECEF')
# plt.show()

