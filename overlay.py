import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import inv


def calculateDifference(pose1, pose2):
    pose1 = np.array(pose1).reshape(4, 4)
    pose2 = np.array(pose2).reshape(4, 4)

    # Extract the rotation matrices and translation vectors from both poses
    R1 = pose1[:3, :3]
    R2 = pose2[:3, :3]
    t1 = pose1[:3, 3]
    t2 = pose2[:3, 3]

    # Calculate the rotation difference between the two poses
    # R_diff = np.linalg.inv(R1) @ R2

    # # Calculate the translation difference between the two poses
    # t_diff = t2 - t1

    # # Print the rotation and translation differences
    # print("Rotation difference:")
    # print(R_diff)

    # print("Translation difference:")
    # print(t_diff)
    # RT = np.hstack((R_diff, t_diff.reshape(3, 1)))
    # RT_homogeneous = np.vstack((RT, np.array([0, 0, 0, 1])))

    # print(RT_homogeneous)
    RT_homogeneous = np.linalg.inv(pose1.T) @ pose2.T
    return RT_homogeneous


f = open('2E8519AA-4E63-4B9E-9D4E-05F7AFA202872023-03-03T09:51:12-05:00-0_metadata.json')
metadata = json.load(f)
f2 = open('2E8519AA-4E63-4B9E-9D4E-05F7AFA202872023-03-03T09:51:12-05:00_pathdata.json')
pathdata = json.load(f2)

alldata = pathdata
for key in metadata:
    alldata[key] = metadata[key]


lat_long = np.array(list(map(lambda x: (x['GARAnchorUUID'], x['longitude'], x['latitude'], x['altitude'], x['geoAnchorTransform']), alldata['savedRouteGeospatialLocations'])))
lat_long_only = np.array(list(map(lambda x: (x[1],x[2],x[3]), lat_long)))
lat_long_poses = np.array(list(map(lambda x: (x[-1]), lat_long)))
gar_anchors = np.array(list(map(lambda x: (x['identifier'], x['transform']), alldata['garAnchors'][-1])))

# filter and sort gar_anchors to match corresponding ID for lat_long points
gar_anchors_filtered = []
for x in lat_long:
  for y in gar_anchors:
    if x[0] == y[0]:  # if IDs match
      gar_anchors_filtered.append(y)
gar_anchors_poses_filtered = np.array(list(map(lambda x: (x[-1]), gar_anchors_filtered)))
df = pd.DataFrame(lat_long, columns=['ID', 'longitude', 'latitude', 'altitude', 'geoAnchorTransform'])
df = df.rename(columns={'geoAnchorTransform': 'pose'}) # Rename the column to match the DataFrame
pose_transforms = np.zeros((len(gar_anchors_poses_filtered), 4, 4))
for i in range(len(gar_anchors_poses_filtered)):
    pose_transforms[i] = calculateDifference(gar_anchors_poses_filtered[i], lat_long_poses[i])
print(pose_transforms)
# Use the first pose transform to transform lat_long_only
transformed_points = np.zeros_like(lat_long_only)
for i in range(len(lat_long_only)):
    homogenous_coord = np.hstack((lat_long_only[i], 1))
    transformed_coord = pose_transforms[1] @ homogenous_coord
    transformed_points[i] = transformed_coord[:3]

# Plot the transformed points
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(transformed_points[:,0], transformed_points[:,1], s=10)

plt.show()