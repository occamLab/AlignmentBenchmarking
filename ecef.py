import math
import json
def llh_to_ecef(lat, lon, alt, hor_unc, alt_unc):
    # Constants
    a = 6378137.0  # Semi-major axis of the Earth (meters)
    f = 1 / 298.257223563  # Flattening factor of the Earth
    e2 = 2 * f - f ** 2  # Square of eccentricity

    # Convert lat and lon from degrees to radians
    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)

    # Calculate the prime vertical radius of curvature (N)
    N = a / math.sqrt(1 - e2 * math.sin(lat_rad) ** 2)

    # Calculate ECEF coordinates
    x = (N + alt) * math.cos(lat_rad) * math.cos(lon_rad)
    y = (N + alt) * math.cos(lat_rad) * math.sin(lon_rad)
    z = ((1 - e2) * N + alt) * math.sin(lat_rad)

    uncertainty_ecef = np.sqrt(hor_unc**2 + alt_unc**2)

    return x, y, z, uncertainty_ecef
f = open('bad_metadata.json')
metadata = json.load(f)
f2 = open('bad_pathdata.json')
pathdata = json.load(f2)

alldata = pathdata
for key in metadata:
    alldata[key] = metadata[key]
    
print(len(alldata["geoSpatialTransformTimes"]))
print(len(alldata["savedRouteGeospatialLocations"]))
print(len(alldata["geoSpatialTransforms"]))
print(len(alldata["garAnchorCameraWorldTransformsAndGeoSpatialData"]))
#print([x - alldata["garAnchorTimestamps"][0] for x in alldata["garAnchorTimestamps"]])
#print([x - alldata["geoSpatialTransformTimes"][0] for x in alldata["geoSpatialTransformTimes"]])