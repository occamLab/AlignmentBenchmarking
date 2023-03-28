from scipy.spatial.transform import Rotation as R

def get_matching_cloud_anchor(cloud_identifier):
    for anchor in alldata['cloudAnchorsForAlignment']:
        if anchor['cloudAnchorID'] == cloud_identifier:
            return np.array(anchor['anchorTransform']).reshape((4,4)).T
    return None

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
        print("flattened", flattened_rotation)
        print("original", alignment_transform)


