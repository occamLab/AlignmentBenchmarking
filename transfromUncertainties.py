import numpy as np
import math
from scipy import stats
# Plot a lot of things
# Sample a bunch of points plug into pca two axis


def make_circle(lat, lon, radius, points,):
    lat = []
    lon = []
    step_size = 2*math.pi/points
    step = 0
    for _ in range(points):
        round(math.cos(step) + radius, 2)
        lat.append(round(math.cos(step) * radius + lat, 2))
        lon.append(round(math.sin(step) * radius + lon, 2))
        step = step + step_size
    return lat, lon


def normal_dist(mean, std, points):
    # look into this more
    # math.randomution \
    # choose random value from distrib
    dist = []
    snd = stats.norm(mean, std)
   # x = np.linspace(mean - std *10, mean + std *10, points)
    dist = snd.pdf(x)
    return dist


def meters_to_lat_long(lat, lon, radius, angle):
    r_earth = 6378100
    dy = radius * math.sin(angle)
    dx = radius * math.cos(angle)
    new_latitude = lat + (dy / r_earth) * (180 / math.pi)
    new_longitude = lon + (dx / r_earth) * \
        (180 / math.pi) / math.cos(lat * math.pi/180)
    return new_latitude, new_longitude


def uniform_dist(start, width, points):
    data_uniform = []
    n = 10000
    data_uniform = stats.uniform.rvs(size=points, loc=start, scale=width)
    return data_uniform


def generate_cloud(lat, lon, alt, heading, ho_uncertainty, he_uncertainty, a_uncertainty, points):
    alt_dist = []
    heading_dist = []
    angles_dist = []
    lat_dist = []
    long_dist = []
    alt_dist = normal_dist(alt, a_uncertainty, points)
    heading_dist = normal_dist(heading, he_uncertainty, points)
    angles_dist = uniform_dist(0, 2*math.pi, points)
    distance_dist = normal_dist([0, 0], ho_uncertainty, points)
    for point in range(points):
        lat_dist, long_dist = meters_to_lat_long(
            lat, lon, distance_dist[point], angles_dist[point])
    return lat_dist, long_dist, alt_dist, heading_dist
    # generate altitude
    # generate heading
    # generate angles
