from scipy.spatial.transform import Rotation as R
import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import inv

def get_matching_cloud_anchor(cloud_identifier):
    for anchor in alldata['cloudAnchorsForAlignment']:
        if anchor['cloudAnchorID'] == cloud_identifier:
            return np.array(anchor['anchorTransform']).reshape((4,4)).T
    return None
f = open('2E8519AA-4E63-4B9E-9D4E-05F7AFA202872023-03-03T09:51:12-05:00-0_metadata.json')
metadata = json.load(f)
f2 = open('2E8519AA-4E63-4B9E-9D4E-05F7AFA202872023-03-03T09:51:12-05:00_pathdata.json')
pathdata = json.load(f2)

alldata = pathdata
for key in metadata:
    alldata[key] = metadata[key]
for anchor in alldata['garAnchors'][-3]:
    if anchor['cloudIdentifier']:
        corresponding_pose = get_matching_cloud_anchor(anchor['cloudIdentifier'])
        if corresponding_pose is None:
            continue
        current_pose = np.array(anchor['transform']).reshape((4,4)).T
        alignment_transform = current_pose @ np.linalg.inv(corresponding_pose)
        v = np.cross(alignment_transform[:-1,1], np.array([0.0, 1.0, 0.0]))
        u = v / np.linalg.norm(v)
        theta = np.arcsin(np.linalg.norm(v))
        flattened_rotation = R.from_rotvec(u*theta).as_matrix() @ alignment_transform[:-1,:-1]
        flattened_transform = np.copy(alignment_transform)
        flattened_transform[:-1, :-1] = flattened_rotation
        print("original", alignment_transform)
transformed_poses = []
lat_long = np.array(list(map(lambda x: (x['GARAnchorUUID'], x['longitude'], x['latitude'], x['altitude'], x['geoAnchorTransform']), alldata['savedRouteGeospatialLocations'])))
lat_long_only = np.array(list(map(lambda x: (x[1],x[2],x[3]), lat_long)))
lat_long_poses = np.array(list(map(lambda x: (x[-1]), lat_long)))
for pose in lat_long_poses:
    # Convert pose to homogeneous transformation matrix
    pose_matrix = np.array(pose).reshape((4, 4)).T
    # Transform pose to new coordinate frame using alignment_transform
    transformed_pose = alignment_transform @ pose_matrix
    # Append transformed pose to list of transformed poses
    transformed_poses.append(transformed_pose)

# Convert list of transformed poses to a 2D numpy array
transformed_poses = np.array(transformed_poses)

# Plot x,y coordinates of transformed poses
plt.plot(transformed_poses[:, 0, 3], transformed_poses[:, 1, 3])
plt.show()
