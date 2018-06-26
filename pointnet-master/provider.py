import os
import sys
import numpy as np
import h5py
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(BASE_DIR)

# Download dataset for point cloud classification
DATA_DIR = os.path.join(BASE_DIR, 'data')
if not os.path.exists(DATA_DIR):
    os.mkdir(DATA_DIR)
if not os.path.exists(os.path.join(DATA_DIR, 'modelnet40_ply_hdf5_2048')):
    www = 'https://shapenet.cs.stanford.edu/media/modelnet40_ply_hdf5_2048.zip'
    zipfile = os.path.basename(www)
    os.system('wget %s; unzip %s' % (www, zipfile))
    os.system('mv %s %s' % (zipfile[:-4], DATA_DIR))
    os.system('rm %s' % (zipfile))


def fix_data_size(data, max_size):
    data_new = data
    for i in range(len(data)):
        len_data = len(data[i])
        if len_data < max_size:
            for j in range(max_size - len_data):
                data_new[i].append([0,0,0])

    return np.array(data_new)

def shuffle_data(data, labels):
    """ Shuffle data and labels.
        Input:
          data: B,N,... numpy array
          label: B,... numpy array
        Return:
          shuffled data, label and shuffle indices
    """
    idx = np.arange(len(labels))
    np.random.shuffle(idx)
    data_new = [data[i] for i in idx]
    labels_new = [labels[i] for i in idx]
    return data_new, labels_new, idx


def rotate_point_cloud(batch_data):
    """ Randomly rotate the point clouds to augument the dataset
        rotation is per shape based along up direction
        Input:
          BxNx3 array, original batch of point clouds
        Return:
          BxNx3 array, rotated batch of point clouds
    """
    
    
    batch_data = np.array(batch_data)
    rotated_data = []
    
    '''
    for k in range(batch_data.shape[0]):
        rotation_angle = np.random.uniform() * 2 * np.pi
        cosval = np.cos(rotation_angle)
        sinval = np.sin(rotation_angle)
        rotation_matrix = np.array([[cosval, 0, sinval],
                                    [0, 1, 0],
                                    [-sinval, 0, cosval]])
        shape_pc = np.array(batch_data[k])
        rotated_data.append(np.dot(shape_pc.reshape((-1, 3)), rotation_matrix).tolist())
    return np.array(rotated_data)
    '''
    
    # GL(3,Z) transformation of the set
    option = np.random.randint(2)
    if option == 0:
        # Rotation
        option = np.random.randint(3)
        if option == 0:
            # Around x-axis
            mat = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
        elif option == 1:
            # Around y-axis
            mat = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
        elif option == 2:
            # Around z-axis
            mat = np.array([[0, 1, 0], [0, 0, -1], [0, 0, 1]])
    elif option == 1:
        # Reflection
        option = np.random.randint(3)
        if option == 0:
            # Around yz-plane
            mat = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
        elif option == 1:
            # Around xz-plane
            mat = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])
        elif option == 2:
            # Around xy-plane
            mat = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
    
    for k in range(batch_data.shape[0]):
        shape_pc = np.array(batch_data[k])
        rotated_data.append(np.dot(shape_pc.reshape((-1, 3)), mat).tolist())
    
    return np.array(rotated_data)


def rotate_point_cloud_by_angle(batch_data, rotation_angle):
    """ Rotate the point cloud along up direction with certain angle.
        Input:
          BxNx3 array, original batch of point clouds
        Return:
          BxNx3 array, rotated batch of point clouds
    """
    rotated_data = np.zeros(batch_data.shape, dtype=np.float32)
    for k in range(batch_data.shape[0]):
        #rotation_angle = np.random.uniform() * 2 * np.pi
        cosval = np.cos(rotation_angle)
        sinval = np.sin(rotation_angle)
        rotation_matrix = np.array([[cosval, 0, sinval],
                                    [0, 1, 0],
                                    [-sinval, 0, cosval]])
        shape_pc = batch_data[k, ...]
        rotated_data[k, ...] = np.dot(shape_pc.reshape((-1, 3)), rotation_matrix)
    return rotated_data


def jitter_point_cloud(batch_data, sigma=0.01, clip=0.05):
    """ Randomly jitter points. jittering is per point.
        Input:
          BxNx3 array, original batch of point clouds
        Return:
          BxNx3 array, jittered batch of point clouds
    """
    B, N, C = batch_data.shape
    assert(clip > 0)
    jittered_data = np.clip(sigma * np.random.randn(B, N, C), -1*clip, clip)
    jittered_data += batch_data
    return jittered_data

def getDataFiles(list_filename):
    return [line.rstrip() for line in open(list_filename)]

def load_h5(h5_filename):
    f = h5py.File(h5_filename)
    data = f['data'][:]
    label = f['label'][:]
    return (data, label)

def loadDataFile(input_list):
    len_input = len(input_list)
    print(len_input)
    data = [input_list[i][0] for i in range(len_input)]
    label = [input_list[i][1] for i in range(len_input)]
    return (data, label)
    #return load_h5(filename)
    
def to_grid(batch_data, size):
    batch_grid = []
    batch_num = batch_data.shape[0]
    num_pts = batch_data.shape[1]
    for i in range(batch_num):
        grid = np.zeros((2*size, 2*size, 2*size))
        data = batch_data[i]
        for k in range(num_pts):
            if isinstance(data[k][0], int) and isinstance(data[k][1], int) and isinstance(data[k][2], int):
                for h in range(3):
                    x = data[k][0] + size
                    y = data[k][1] + size
                    z = data[k][2] + size
                    grid[x][y][z] = 1
        grid = grid.flatten()
        batch_grid.append(grid)
    batch_grid = np.array(batch_grid)
    return batch_grid
    
def load_h5_data_label_seg(h5_filename):
    f = h5py.File(h5_filename)
    data = f['data'][:]
    label = f['label'][:]
    seg = f['pid'][:]
    return (data, label, seg)


def loadDataFile_with_seg(filename):
    return load_h5_data_label_seg(filename)
