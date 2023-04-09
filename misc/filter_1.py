import numpy as np
import math
from filterpy.monte_carlo import systematic_resample
import json
class Particle:
    def __init__(self, lat, lon, alt, heading, weight=1.0):
        self.lat = lat
        self.lon = lon
        self.alt = alt
        self.heading = heading
        self.weight = weight

def estimate_new_position(lat, lon, heading, distance):
    R = 6371.0  # Earth radius in km

    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)
    heading_rad = math.radians(heading)

    angular_distance = distance / R

    new_lat_rad = math.asin(
        math.sin(lat_rad) * math.cos(angular_distance)
        + math.cos(lat_rad) * math.sin(angular_distance) * math.cos(heading_rad)
    )
    new_lat = math.degrees(new_lat_rad)

    new_lon_rad = lon_rad + math.atan2(
        math.sin(heading_rad) * math.sin(angular_distance) * math.cos(lat_rad),
        math.cos(angular_distance) - math.sin(lat_rad) * math.sin(new_lat_rad),
    )
    new_lon = math.degrees(new_lon_rad)

    return new_lat, new_lon

def predict(particles, delta_time, std_devs, velocity):
    for p in particles:
        distance_traveled = velocity * delta_time
        new_lat, new_lon = estimate_new_position(p.lat, p.lon, p.heading, distance_traveled)
        p.lat = new_lat + std_devs[0] * np.random.randn()
        p.lon = new_lon + std_devs[1] * np.random.randn()
        p.alt += delta_time * std_devs[2] * np.random.randn()
        p.heading += delta_time * std_devs[3] * np.random.randn()

def update(particles, garAnchors, measurement, std_devs):
    for p in particles:
        if measurement["GARAnchorUUID"] != "":
            anchor = garAnchors[measurement["GARAnchorUUID"]]
            p.weight *= np.exp(
                -0.5 * (
                    ((p.lat - anchor["latitude"]) / std_devs[0]) ** 2
                    + ((p.lon - anchor["longitude"]) / std_devs[1]) ** 2
                    + ((p.alt - anchor["altitude"]) / std_devs[2]) ** 2
                    + ((p.heading - anchor["heading"]) / std_devs[3]) ** 2
                )
            )
        else:
            p.weight *= np.exp(
                -0.5 * (
                    ((p.lat - measurement["latitude"]) / std_devs[0]) ** 2
                    + ((p.lon - measurement["longitude"]) / std_devs[1]) ** 2
                    + ((p.alt - measurement["altitude"]) / std_devs[2]) ** 2
                    + ((p.heading - measurement["heading"]) / std_devs[3]) ** 2
                )
            )
    normalize_weights(particles)

def normalize_weights(particles):
    weights = np.array([p.weight for p in particles])
    eps = 1e-10
    weights /= (np.sum(weights) + eps)
    for i, p in enumerate(particles):
        p.weight = weights[i]


def multinomial_resample(weights):
    N = len(weights)
    indices = np.zeros(N, dtype=int)
    cdf = np.cumsum(weights)
    
    for i in range(N):
        random_sample = np.random.random()
        index = np.searchsorted(cdf, random_sample)
        indices[i] = index
        
    return indices


def resample(particles):
    weights = np.array([p.weight for p in particles])
    indices = multinomial_resample(weights)
    print(f"indices: {indices}")
    print(f"max index: {max(indices)}, particles length: {len(particles)}")
    return resample_from_index(particles, indices)


def resample_from_index(particles, indices):
    try:
        new_particles = [particles[i] for i in indices]
    except IndexError as e:
        print(f"Error in resampling: {e}")
        print(f"indices: {indices}")
        print(f"particles length: {len(particles)}")
        raise e

    for p in new_particles:
        p.weight = 1 / len(new_particles)
    return new_particles





def estimate(particles):
    lat = np.mean([p.lat for p in particles])
    lon = np.mean([p.lon for p in particles])
    alt = np.mean([p.alt for p in particles])
    heading = np.mean([p.heading for p in particles])
    return lat, lon, alt, heading

def particle_filter(savedRouteGeoSpatialTransforms, garAnchors, num_particles=1000, velocity=0.1):
    particles = [ Particle(np.random.normal(savedRouteGeoSpatialTransforms[0]["latitude"], savedRouteGeoSpatialTransforms[0]["horizontalUncertainty"]),
    np.random.normal(savedRouteGeoSpatialTransforms[0]["longitude"], savedRouteGeoSpatialTransforms[0]["horizontalUncertainty"]),
    np.random.normal(savedRouteGeoSpatialTransforms[0]["altitude"], savedRouteGeoSpatialTransforms[0]["altitudeUncertainty"]),
    np.random.normal(savedRouteGeoSpatialTransforms[0]["heading"], savedRouteGeoSpatialTransforms[0]["headingUncertainty"]),
    )
    for _ in range(num_particles)
    ]
    process_std_devs = (0.1, 0.1, 0.1, 1.0)
    measurement_std_devs = (
        savedRouteGeoSpatialTransforms[0]["horizontalUncertainty"],
        savedRouteGeoSpatialTransforms[0]["horizontalUncertainty"],
        savedRouteGeoSpatialTransforms[0]["altitudeUncertainty"],
        savedRouteGeoSpatialTransforms[0]["headingUncertainty"],
    )

    estimated_states = []

    for i in range(1, len(savedRouteGeoSpatialTransforms)):
        delta_time = 0.1

        predict(particles, delta_time, process_std_devs, velocity)
        update(particles, garAnchors, savedRouteGeoSpatialTransforms[i], measurement_std_devs)

        particles = resample(particles)
        estimated_state = estimate(particles)
        estimated_states.append(estimated_state)

    return estimated_states
f = open('bad_metadata.json')
metadata = json.load(f)
f2 = open('bad_pathdata.json')
pathdata = json.load(f2)

alldata = pathdata
for key in metadata:
    alldata[key] = metadata[key]
savedRouteGeoSpatialTransforms = alldata["savedRouteGeospatialLocations"]

geoLocationAlignmentAttempts = alldata["geoLocationAlignmentAttempts"]

garAnchors = {
    attempt["geoSpatialAlignmentCrumb"]["GARAnchorUUID"]: {
        "latitude": attempt["geoSpatialAlignmentCrumb"]["latitude"],
        "longitude": attempt["geoSpatialAlignmentCrumb"]["longitude"],
        "altitude": attempt["geoSpatialAlignmentCrumb"]["altitude"],
        "heading": attempt["geoSpatialAlignmentCrumb"]["heading"],
        "horizontalUncertainty": attempt["geoSpatialAlignmentCrumb"]["horizontalUncertainty"],
        "headingUncertainty": attempt["geoSpatialAlignmentCrumb"]["headingUncertainty"],
        "altitudeUncertainty": attempt["geoSpatialAlignmentCrumb"]["altitudeUncertainty"],
    }
    for attempt in geoLocationAlignmentAttempts
    if attempt["wasAccepted"]
}

estimated_states = particle_filter(savedRouteGeoSpatialTransforms, garAnchors)

print("Estimated states:")
for state in estimated_states:
    print(f"Latitude: {state[0]}, Longitude: {state[1]}, Altitude: {state[2]}, Heading: {state[3]}")


#print(len(alldata["geoSpatialTransformTimes"]))

#print(len(alldata["geoSpatialTransforms"]))
#print(len(alldata["garAnchorCameraWorldTransformsAndGeoSpatialData"]))
#print([x - alldata["garAnchorTimestamps"][0] for x in alldata["garAnchorTimestamps"]])
#print([x - alldata["geoSpatialTransformTimes"][0] for x in alldata["geoSpatialTransformTimes"]])
