import math

def llh_to_ecef(lat, lon, alt):
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

    return x, y, z
