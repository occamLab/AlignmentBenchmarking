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
from scipy.spatial.transform import Rotation as R


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
s

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
    pose = np.array(d["cameraWorldTransform"]).reshape(4, 4).T
    esu = d["geospatialTransform"]['eastUpSounth']
    coords.append([lat, lon, alt, heading, lat_uncertainty, lon_uncertainty, alt_uncertainty, pose, esu])

transform_poses = [np.array(x["cameraWorldTransform"]).reshape(4, 4).T for x in alldata["garAnchorCameraWorldTransformsAndGeoSpatialData"]]
inv(transform_poses[1]) @ transform_poses[0]

predictions = particle_filter(coords)

ecef_coords = [latlonalt_to_ecef(coord[0], coord[1], coord[2]) for coord in coords]
x_coords = [coord[0] for coord in ecef_coords]
y_coords = [coord[1] for coord in ecef_coords]

plt.plot(x_coords, y_coords, 'go', label='Ground Truth')
plt.plot([pred[0] for pred in predictions], [pred[1] for pred in predictions], 'bx', label='Predicted')
plt.legend(loc='lower right')
plt.xlabel('ECEF X')
plt.ylabel('ECEF Y')
plt.title('Particle Filter Predictions in ECEF')
plt.show()

