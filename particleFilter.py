import numpy as np
from filterpy.monte_carlo import systematic_resample
from numpy.linalg import norm
from numpy.random import randn
import matplotlib.pyplot as plt
import json

def create_particles(geoSpatialTransforms, N):
    particles = np.empty((N, 3))
    lat, lon, alt = geoSpatialTransforms[0][:3]
    particles[:, 0] = lat
    particles[:, 1] = lon
    particles[:, 2] = alt
    print(particles)
    return particles

def predict(particles, delta_time, heading, speed=0.0001):
    lat, lon, _ = particles.T
    delta_lat = np.sin(heading) * speed * delta_time
    delta_lon = np.cos(heading) * speed * delta_time
    particles[:, 0] += delta_lat
    particles[:, 1] += delta_lon

# def predict(particles, delta_time, heading, velocity):
#     for p in particles:
#         distance_traveled = velocity * delta_time
#         new_lat, new_lon = estimate_new_position(p.lat, p.lon, p.heading, distance_traveled)
#         p.lat = new_lat + std_devs[0] * np.random.randn()
#         p.lon = new_lon + std_devs[1] * np.random.randn()
#         p.alt += delta_time * std_devs[2] * np.random.randn()
#         p.heading += delta_time * std_devs[3] * np.random.randn()
def update(particles, weights, measurement, uncertainties):
    lat_uncertainty, lon_uncertainty, alt_uncertainty = uncertainties
    weights.fill(1.0)

    for i, p in enumerate(particles):
        dist = np.sqrt(((p[0] - measurement[0]) / (lat_uncertainty*0.1)) ** 2 +
                       ((p[1] - measurement[1]) / (lon_uncertainty*0.1)) ** 2 +
                       ((p[2] - measurement[2]) / (alt_uncertainty*0.1)) ** 2)
        weights[i] *= np.exp(-dist)

    weights += 1.e-300
    weights /= np.sum(weights)

def estimate(particles, weights):
    return np.average(particles, axis=0, weights=weights)

def neff(weights):
    return 1. / np.sum(np.square(weights))

def resample_from_index(particles, weights, indexes):
    particles[:] = particles[indexes]
    weights.resize(len(particles))
    weights.fill(1.0 / len(weights))

def particle_filter(geoSpatialTransforms, timestamps, N=1000):
    particles = create_particles(geoSpatialTransforms, N)
    weights = np.ones(N) / N

    predictions = []

    for t in range(1, len(geoSpatialTransforms)):
        delta_time = (timestamps[t] - timestamps[t - 1])
        heading = geoSpatialTransforms[t][3]
        predict(particles, delta_time, heading)

        measurement = geoSpatialTransforms[t][:3]
        uncertainties = geoSpatialTransforms[t][4:]
        update(particles, weights, measurement, uncertainties)

        if neff(weights) < N / 2:
            indexes = systematic_resample(weights)
            resample_from_index(particles, weights, indexes)

        mu = estimate(particles, weights)
        predictions.append(mu)

    return predictions

# Example usage
f = open('bad_metadata.json')
metadata = json.load(f)
f2 = open('bad_pathdata.json')
pathdata = json.load(f2)

alldata = pathdata
for key in metadata:
    alldata[key] = metadata[key]
timestamps = alldata["geoSpatialTransformTimes"]
coords = []
for d in alldata['geoSpatialTransforms']:
    lat = d['latitude']
    lon = d['longitude']
    alt = d['altitude']
    heading = d['heading']
    lat_uncertainty = d['positionAccuracy']
    lon_uncertainty = d['positionAccuracy']
    alt_uncertainty = d['altitudeAccuracy']
    coords.append([lat, lon, alt, heading, lat_uncertainty, lon_uncertainty, alt_uncertainty])

predictions = particle_filter(coords, timestamps)

latitudes = [x[0] for x in coords]
longitudes = [x[1] for x in coords]
plt.plot(longitudes, latitudes, 'go', label='Ground Truth')
plt.plot([x[0] for x in predictions], [x[1] for x in predictions], 'bo', label='Predicted')
plt.legend(loc='lower right')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Particle Filter Predictions')
plt.show()
