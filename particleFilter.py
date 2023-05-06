import numpy as np
from filterpy.monte_carlo import systematic_resample
from numpy.linalg import norm
from numpy.random import randn
import matplotlib.pyplot as plt
import json
from pyproj import Transformer
from scipy.stats import multivariate_normal
from scipy.linalg import inv
from scipy.spatial.transform import Rotation as R
from particle_cloud import generate_cloud

def enu_to_ecef_rotation_matrix(lat, lon):
    sin_lat = np.sin(np.radians(lat))
    cos_lat = np.cos(np.radians(lat))
    sin_lon = np.sin(np.radians(lon))
    cos_lon = np.cos(np.radians(lon))

    R = np.array([[-sin_lon, -sin_lat * cos_lon, cos_lat * cos_lon],
                  [cos_lon, -sin_lat * sin_lon, cos_lat * sin_lon],
                  [0, cos_lat, sin_lat]])
    return R

def esu_to_ecef_rotation_matrix(esu):
    axis = np.array(esu['axis'])
    angle = esu['angle']
    return R.from_rotvec(axis*angle).as_matrix()

def create_particles(geoSpatialTransforms, N):
    particles = [dict() for _ in range(N)]
    high_acc_pt = geoSpatialTransforms[0]
    lat, lon, alt, _ = generate_cloud(high_acc_pt[0], high_acc_pt[1], high_acc_pt[2], high_acc_pt[3], high_acc_pt[4], high_acc_pt[5], high_acc_pt[6], N)
    x_coords, y_coords, z_coords = [], [], []
    for i in range(N):
        x, y, z = latlonalt_to_ecef(lat[i], lon[i], alt[i])
        x_coords.append(x)
        y_coords.append(y)
        z_coords.append(z)
    pose = geoSpatialTransforms[0][-2]
    for i, particle in enumerate(particles):
        particle['x'] = x_coords[i]
        particle['y'] = y_coords[i]
        particle['z'] = z_coords[i]
        particle['pose'] = pose
    return particles

def latlonalt_to_ecef(lat, lon, alt):
    transformer = Transformer.from_crs("epsg:4326", "epsg:4978")
    x, y, z = transformer.transform(lat, lon, alt)
    return x, y, z

def ar_pose_to_ecef(ar_pose, esu, lat, lon, alt):
    ar_pose_matrix = np.array(ar_pose).reshape(4, 4)

    # Convert LLA to ECEF
    x, y, z = latlonalt_to_ecef(lat, lon, alt)

    # Get the rotation matrix from ESU to ECEF
    R = esu_to_ecef_rotation_matrix(esu)

    # Construct the ECEF transformation matrix
    T_ecef = np.eye(4)
    T_ecef[:3, :3] = R
    T_ecef[:3, 3] = [x, y, z]

    # Convert AR pose from ESU to ECEF
    ar_pose_ecef = np.dot(T_ecef, ar_pose_matrix)

    return ar_pose_ecef

def update_ecef_uncertainties(lat, lon, alt, lat_uncertainty, lon_uncertainty, alt_uncertainty):
    a = 6378137.0  # Earth's semi-major axis (m)
    f = 1 / 298.257223563  # Earth's flattening
    e2 = 2 * f - f ** 2  # Earth's first eccentricity squared
    # Not very good
    N_phi = a / np.sqrt(1 - e2 * np.sin(np.radians(lat)) ** 2)
    x_uncertainty = np.abs((N_phi + alt) * np.cos(np.radians(lat)) * np.radians(lon_uncertainty))
    y_uncertainty = np.abs((N_phi * (1 - e2) + alt) * np.sin(np.radians(lat)) * np.radians(lat_uncertainty))
    z_uncertainty = alt_uncertainty

    return x_uncertainty, y_uncertainty, z_uncertainty
def predict(particles, prev_lat, prev_lon, prev_alt, current_lat, current_lon, current_alt, prev_pose, pose, prev_esu, esu, noise_std=1):
    prev_x, prev_y, prev_z = latlonalt_to_ecef(prev_lat, prev_lon, prev_alt)
    current_x, current_y, current_z = latlonalt_to_ecef(current_lat, current_lon, current_alt)
    dist_ecef = np.array([current_x - prev_x, current_y - prev_y, current_z - prev_z])

    pose = ar_pose_to_ecef(pose, esu, current_lat, current_lon, current_alt)
    prev_pose = ar_pose_to_ecef(prev_pose, prev_esu, prev_lat, prev_lon, prev_alt)

    pose_difference = inv(prev_pose) @ pose

    translation = pose_difference[0:3, 3]
    # rotation = pose_difference[0:3, 0:3]
    noise = np.random.normal(0, noise_std, size=(len(particles), 3))
    # for i, particle in enumerate(particles):
    #     particle['x'] += dist_ecef[0] + noise[i, 0]
    #     particle['y'] += dist_ecef[1] + noise[i, 1]
    #     particle['z'] += dist_ecef[2] + noise[i, 2]
    for i, particle in enumerate(particles):
        particle['x'] += translation[0]
        particle['y'] += translation[1]
        particle['z'] += translation[2]
        particle['pose'] = inv(pose_difference) @ particle['pose'] 

def update(particles, weights, measurement_ecef, uncertainties, weight_effect=1):
    weights.fill(1.0)

    for i, p in enumerate(particles):
        particle_array = np.array([p['x'], p['y'], p['z']])
        cov = np.diag(np.array(uncertainties) ** 2)
        weight = multivariate_normal.pdf(particle_array, mean=measurement_ecef, cov=cov)
        weights[i] = weight ** weight_effect
    weights += 1.e-300
    weights /= np.sum(weights)
    return weights

def estimate(particles, weights):
    x = np.average([p['x'] for p in particles], weights=weights)
    y = np.average([p['y'] for p in particles], weights=weights)
    z = np.average([p['z'] for p in particles], weights=weights)
    return np.array([x, y, z])

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
        prev_pose = geoSpatialTransforms[t-1][-2]
        pose = geoSpatialTransforms[t][-2]
        prev_esu = geoSpatialTransforms[t-1][-1]
        esu = geoSpatialTransforms[t][-1]
        current_lat, current_lon, current_alt = geoSpatialTransforms[t][:3]
        heading = geoSpatialTransforms[t][3]
        #predict(particles, delta_time, prev_lat, prev_lon, prev_alt, heading)
        predict(particles, prev_lat, prev_lon, prev_alt, current_lat, current_lon, current_alt, prev_pose, pose, prev_esu, esu)
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
    alt_uncertainty = d["geospatialTransform"]['altitudeAccuracy']
    heading_uncertainty = d["geospatialTransform"]['orientationYawAccuracy']
    pose = np.array(d["cameraWorldTransform"]).reshape(4, 4).T
    esu = d["geospatialTransform"]['eastUpSounth']
    coords.append([lat, lon, alt, heading, lat_uncertainty, alt_uncertainty, heading_uncertainty, pose, esu])

transform_poses = [np.array(x["cameraWorldTransform"]).reshape(4, 4).T for x in alldata["garAnchorCameraWorldTransformsAndGeoSpatialData"]]
inv(transform_poses[1]) @ transform_poses[0]

predictions = particle_filter(coords)

ecef_coords = [latlonalt_to_ecef(coord[0], coord[1], coord[2]) for coord in coords]
x_coords = [coord[0] for coord in ecef_coords]
y_coords = [coord[1] for coord in ecef_coords]
z_coords = [coord[2] for coord in ecef_coords]

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter(x_coords, y_coords, z_coords, c='g', marker='o', label='Ground Truth')
ax.scatter([pred[0] for pred in predictions], [pred[1] for pred in predictions], [pred[2] for pred in predictions], c='b', marker='x', label='Predicted')
# plt.plot(x_coords, y_coords, z_coords, 'go', label='Ground Truth')
# plt.plot([pred[0] for pred in predictions], [pred[1] for pred in predictions], [pred[2] for pred in predictions],'bx', label='Predicted')

# plt.legend(loc='lower right')
# plt.xlabel('ECEF X')
# plt.ylabel('ECEF Y')
plt.title('Particle Filter Predictions in ECEF')
plt.show()

