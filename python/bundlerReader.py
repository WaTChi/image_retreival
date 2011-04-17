def read_bundle_output(filename):
    with open(filename) as file:
        # first line is a useless comment
        file.readline()
        num_cameras, num_points = [int(num) for num in file.readline().split()]
        cameras = read_camera_info(file, num_cameras)
        points = read_point_info(file, num_points)
        return cameras, points

def read_camera_info(file, num_cameras):
    result = []
    for camera_num in range(num_cameras):
        line = [float(num) for num in file.readline().split()]
        focal_length, radial_dist = line[0], line[1:]
        rotation = [[float(num) for num in file.readline().split()] for line in range(3)]
        translation = [float(num) for num in file.readline().split()]
        result.append((focal_length, radial_dist, rotation, translation))
    return result

def read_point_info(file, num_points):
    result = []
    for point_num in range(num_points):
        position = [float(num) for num in file.readline().split()]
        color = [int(num) for num in file.readline().split()]
        views = [float(num) for num in file.readline().split()]
        result.append((position, color, views))
    return result

