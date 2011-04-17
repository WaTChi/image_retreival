import os.path
import os
import shutil
import bundlerReader
import math

TMPDIR = '/tmp/bundler'
BUNDLERPATH = '/home/etzeng/bundler-v0.3-binary/RunBundler.sh'


def run_bundler(files):
    """
    Run bundler on the list of files given as argument. Leaves a file named
    bundle.out in the current working directory.
    """
    cwd = os.getcwd()
    if os.path.exists(TMPDIR):
        shutil.rmtree(TMPDIR)
    os.mkdir(TMPDIR)
    for filename in files:
        shutil.copy(filename, TMPDIR)
    os.chdir(TMPDIR)
    os.system(BUNDLERPATH)
    outfile = os.path.join(TMPDIR, "bundle/bundle.out")
    if os.path.exists(outfile):
        shutil.copy(outfile, cwd)
        os.chdir(cwd)

def find_r_and_t(file1, file2):
    """
    Uses bundler to find the rotation and translation vectors to take one image
    to the next.
    """
    run_bundler([file1, file2])
    if not os.path.exists("bundle.out"):
        return None, None
    cameras, points = bundlerReader.read_bundle_output("bundle.out")
    r = vector_diff(r_to_euler_angles(cameras[1][2]),
                    r_to_euler_angles(cameras[0][2]))
    t = vector_diff(cameras[1][3], cameras[0][3])
    return [val * 180 / math.pi for val in r], t

def vector_diff(v1, v2):
    """ Perform a 'dot subtract' on two vectors. """
    return [x - y for x, y in zip(v1, v2)]

def r_to_euler_angles(R):
    """ Converts a rotation matrix to pitch, yaw, roll """
    pitch = math.atan2(R[2][0], R[2][1])
    yaw = math.acos(R[2][2])
    roll = -math.atan2(R[0][2], R[1][2])
    return pitch, yaw, roll

