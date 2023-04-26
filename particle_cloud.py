import numpy as np
from numpy.linalg import norm
import math
from numpy.random import randn
import matplotlib.pyplot as plt
import json
from pyproj import Transformer
from scipy.linalg import inv
from scipy.spatial.transform import Rotation as R
from scipy import stats


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
    normal = np.random.normal(mean, std, points)
    dist = np.random.choice(normal, size=points)
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
    distance_dist = normal_dist(0, ho_uncertainty, points)
    for point in range(points):
        lat_dist, long_dist = meters_to_lat_long(
            lat, lon, distance_dist[point], angles_dist[point])
    return lat_dist, long_dist, alt_dist, heading_dist


def create_particles(geoSpatialTransforms, N):
    particles = [dict() for _ in range(N)]
    lat, lon, alt = geoSpatialTransforms[0][:3]
    x, y, z = latlonalt_to_ecef(lat, lon, alt)
    pose = geoSpatialTransforms[0][-2]
    for particle in particles:
        particle['x'] = x
        particle['y'] = y
        particle['z'] = z
        particle['pose'] = pose
    print(particle)
    return particles


def latlonalt_to_ecef(lat, lon, alt):
    transformer = Transformer.from_crs("epsg:4326", "epsg:4978")
    x, y, z = transformer.transform(lat, lon, alt)
    return x, y, z


f = open('eastupsouth_metadata.json')
metadata = json.load(f)
f2 = open('eastupsouth_pathdata.json')
pathdata = json.load(f2)

alldata = pathdatagit
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
    esu = d["geospatialTransform"]['eastUpSounth']
    coords.append([lat, lon, alt, heading, lat_uncertainty,
                  lon_uncertainty, alt_uncertainty, pose, esu])

lat_dist, long_dist, alt_dist, heading_dist = generate_cloud(
    lat, lon, alt, heading, lat_uncertainty, alt_uncertainty, .01, 1000)
