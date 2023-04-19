import numpy as np
import math

def make_circle(origin, radius, points,):
    lat = []
    lon = []
    step_size = 2*math.pi/points
    step = 0
    for _ in range(points):
        round(math.cos(step) + radius, 2)
        lat.append(round(math.cos(step) * radius + origin[1], 2))
        lon.append(round(math.sin(step) * radius + origin[2], 2))
        step = step + step_size
    return lat, lon
