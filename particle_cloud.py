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
import open3d as o3d
from sklearn.decomposition import PCA
import pandas as pd


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


def generate_cloud(lat, lon, alt, heading, ho_uncertainty, a_uncertainty, he_uncertainty, points):
    alt_dist = []
    heading_dist = []
    angles_dist = []
    lat_dist = np.zeros(points)
    long_dist = np.zeros(points)
    alt_dist = normal_dist(alt, a_uncertainty, points)
    heading_dist = normal_dist(heading, he_uncertainty, points)
    angles_dist = uniform_dist(0, 2*math.pi, points)
    distance_dist = normal_dist(0, ho_uncertainty, points)
    for point in range(points):
        lat_dist[point], long_dist[point] = meters_to_lat_long(
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
    yaw_uncertainty = d["geospatialTransform"]['orientationYawAccuracy']
    pose = np.array(d["cameraWorldTransform"]).reshape(4, 4).T
    esu = d["geospatialTransform"]['eastUpSounth']
    coords.append([lat, lon, alt, heading, lat_uncertainty,
                  lon_uncertainty, alt_uncertainty, pose, esu])
points = 1000
print(lat)
print(lon)
print(lat_uncertainty)
print(alt)
print(alt_uncertainty)
print(heading)
print(yaw_uncertainty)
lat_dist, long_dist, alt_dist, heading_dist = generate_cloud(
    lat, lon, alt, heading, lat_uncertainty, alt_uncertainty, yaw_uncertainty, points)

plt.hist(lat_dist, color='lightgreen', ec='black', bins=15)
plt.title("Latitude Distribution")
plt.show()
plt.hist(long_dist, color='lightgreen', ec='black')
plt.title("Longitude Distribution")
plt.show()
plt.hist(heading_dist, color='lightgreen', ec='black')
plt.title("Heading Distribution")
plt.show()
plt.hist(alt_dist, color='lightgreen', ec='black')
plt.title("Altitude Distribution")
plt.show()
plt.show()
plt.hist2d(lat_dist, long_dist)
plt.title("Lat/Long Distribution")
plt.show()
x_dist = np.zeros(points)
y_dist = np.zeros(points)
z_dist = np.zeros(points)
for point in range(points):
    x_dist[point], y_dist[point], z_dist[point] = latlonalt_to_ecef(
        lat_dist[point], long_dist[point], alt_dist[point])


fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter(x_dist, y_dist, z_dist, c=z_dist, cmap='viridis', linewidth=0.5)
plt.title("ECEF Distribution")
plt.show()

ax = plt.axes(projection='3d')
ax.plot_trisurf(x_dist, y_dist, z_dist,
                cmap='viridis', edgecolor='none')
plt.title("ECEF Distribution")
plt.show()

pca = PCA(n_components=2)
lat_lon_dataset = pd.DataFrame({'lat': lat_dist, 'lon': long_dist})
print(lat_lon_dataset)
pca.fit(lat_lon_dataset)
principalComponents = pca.fit_transform(lat_lon_dataset)
principalDf = pd.DataFrame(data=principalComponents, columns=[
                           'principal component 1', 'principal component 2'])
plt.figure(figsize=(10, 10))
plt.xticks(fontsize=12)
plt.yticks(fontsize=14)
plt.xlabel('Principal Component - 1', fontsize=20)
plt.ylabel('Principal Component - 2', fontsize=20)
plt.title("Principal Component Analysis of lat/lon Dataset", fontsize=20)
plt.scatter(principalDf['principal component 1'],
            principalDf['principal component 2'], c='r', s=50)
plt.show()
