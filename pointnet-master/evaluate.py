import tensorflow as tf
import math
import numpy as np
import argparse
import socket
import importlib
import time
import os
import scipy.misc
import sys
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(BASE_DIR)
sys.path.append(os.path.join(BASE_DIR, 'models'))
sys.path.append(os.path.join(BASE_DIR, 'utils'))
import provider
import pc_util


parser = argparse.ArgumentParser()
parser.add_argument('--gpu', type=int, default=0, help='GPU to use [default: GPU 0]')
parser.add_argument('--model', default='pointnet_cls', help='Model name: pointnet_cls or pointnet_cls_basic [default: pointnet_cls]')
parser.add_argument('--batch_size', type=int, default=4, help='Batch Size during training [default: 1]')
parser.add_argument('--num_point', type=int, default=128, help='Point Number [256/512/1024/2048] [default: 128]')
parser.add_argument('--model_path', default='log/model.ckpt', help='model checkpoint file path [default: log/model.ckpt]')
parser.add_argument('--dump_dir', default='dump', help='dump folder path [dump]')
parser.add_argument('--visu', action='store_true', help='Whether to dump image for error case [default: False]')
FLAGS = parser.parse_args()


BATCH_SIZE = FLAGS.batch_size
NUM_POINT = FLAGS.num_point
MODEL_PATH = FLAGS.model_path
GPU_INDEX = FLAGS.gpu
MODEL = importlib.import_module(FLAGS.model) # import network module
DUMP_DIR = FLAGS.dump_dir
if not os.path.exists(DUMP_DIR): os.mkdir(DUMP_DIR)
LOG_FOUT = open(os.path.join(DUMP_DIR, 'log_evaluate.txt'), 'w')
LOG_FOUT.write(str(FLAGS)+'\n')

NUM_CLASSES = 1
SHAPE_NAMES = [line.rstrip() for line in \
    open(os.path.join(BASE_DIR, 'data/modelnet40_ply_hdf5_2048/shape_names.txt'))] 

HOSTNAME = socket.gethostname()

# ModelNet40 official train/test split
#TRAIN_FILES = provider.getDataFiles( \
#    os.path.join(BASE_DIR, 'data/modelnet40_ply_hdf5_2048/train_files.txt'))
TRAIN_FILES_path = '/home/carnd/CYML/output/train/train_cube_5_35.txt'
TRAIN_FILES_path = str(TRAIN_FILES_path)

train_file = open(TRAIN_FILES_path, 'r')

TRAIN_FILES = [ eval(line) for line in train_file]
TRAIN_FILES = np.array(TRAIN_FILES)
print("Train file shape: ", TRAIN_FILES.shape)
TRAIN_FILES = TRAIN_FILES.reshape((6, -1, 2))
print("Train file loaded: ", TRAIN_FILES.size)
TRAIN_FILES = TRAIN_FILES.tolist()
print("Sample: ", TRAIN_FILES[0][0])

#TEST_FILES = provider.getDataFiles(\
#    os.path.join(BASE_DIR, 'data/modelnet40_ply_hdf5_2048/test_files.txt'))

TEST_FILES_path = '/home/carnd/CYML/output/train/cube_3_35.txt'
TEST_FILES_path = str(TEST_FILES_path)

test_file = open(TEST_FILES_path, 'r')

TEST_FILES = [eval(line) for line in test_file]
TEST_FILES = np.array(TEST_FILES)
TEST_FILES = TEST_FILES.reshape((2, -1, 2))
print("Test file loaded: ", TEST_FILES.size)
TEST_FILES = TEST_FILES.tolist()
print("Sample: ", TEST_FILES[0][0])

def log_string(out_str):
    LOG_FOUT.write(out_str+'\n')
    LOG_FOUT.flush()
    print(out_str)

def evaluate(num_votes):
    is_training = False
     
    with tf.device('/gpu:'+str(GPU_INDEX)):
        pointclouds_pl, labels_pl = MODEL.placeholder_inputs(BATCH_SIZE, NUM_POINT)
        is_training_pl = tf.placeholder(tf.bool, shape=())

        # simple model
        pred, end_points = MODEL.get_model(pointclouds_pl, is_training_pl)
        loss = MODEL.get_loss(pred, labels_pl, end_points)
        
        # Add ops to save and restore all the variables.
        saver = tf.train.Saver()
        
    # Create a session
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    config.allow_soft_placement = True
    config.log_device_placement = True
    sess = tf.Session(config=config)

    # Restore variables from disk.
    saver.restore(sess, MODEL_PATH)
    log_string("Model restored.")

    ops = {'pointclouds_pl': pointclouds_pl,
           'labels_pl': labels_pl,
           'is_training_pl': is_training_pl,
           'pred': pred,
           'loss': loss}
    
    out_path = 'dump/error_pts_6_25_1545.txt'
    
    print("Sample Test: ", TEST_FILES[0][0])
    
    eval_one_epoch(sess, ops, out_path, num_votes)

   
def eval_one_epoch(sess, ops, out_path, num_votes=1):
    error_cnt = 0
    is_training = False
    total_correct = 0
    total_seen = 0
    loss_sum = 0
    #total_seen_class = [0 for _ in range(NUM_CLASSES)]
    #total_correct_class = [0 for _ in range(NUM_CLASSES)]
    fout = open(os.path.join(DUMP_DIR, 'pred_label.txt'), 'w')
    for fn in range(len(TEST_FILES)):
        log_string('----'+str(fn)+'----')
        current_data, current_label = provider.loadDataFile(TEST_FILES[fn])
        current_data = np.array(current_data)
        current_label = np.array(current_label)
        current_label = np.squeeze(current_label)
        print(current_data.shape)
        
        file_size = current_data.shape[0]
        num_batches = math.floor(file_size // BATCH_SIZE)
        print(file_size)
        
        for batch_idx in range(num_batches):
            start_idx = batch_idx * BATCH_SIZE
            end_idx = (batch_idx+1) * BATCH_SIZE
            cur_batch_size = end_idx - start_idx
            
            # Aggregating BEG
            batch_loss_sum = 0 # sum of losses for the batch
            #batch_pred_sum = np.zeros((cur_batch_size, NUM_CLASSES)) # score for classes
            #batch_pred_classes = np.zeros((cur_batch_size, NUM_CLASSES)) # 0/1 for classes
            for vote_idx in range(num_votes):
                rotated_data = provider.rotate_point_cloud_by_angle(current_data[start_idx:end_idx],
                                                  vote_idx/float(num_votes) * np.pi * 2)
                feed_dict = {ops['pointclouds_pl']: rotated_data,
                             ops['labels_pl']: current_label[start_idx:end_idx],
                             ops['is_training_pl']: is_training}
                loss_val, pred_val = sess.run([ops['loss'], ops['pred']],
                                          feed_dict=feed_dict)
                
                print("input_data: ", rotated_data)
                print("input_label: ", current_label[start_idx:end_idx])
                print("pred_val: ", pred_val)
                #batch_pred_sum += pred_val
                #batch_pred_val = np.argmax(pred_val, 1)
                #for el_idx in range(cur_batch_size):
                #    batch_pred_classes[el_idx, batch_pred_val[el_idx]] += 1
                batch_loss_sum += (loss_val * cur_batch_size / float(num_votes))
            # pred_val_topk = np.argsort(batch_pred_sum, axis=-1)[:,-1*np.array(range(topk))-1]
            # pred_val = np.argmax(batch_pred_classes, 1)
            #pred_val = np.argmax(batch_pred_sum, 1) #return the index of the highest prediction
            # Aggregating END
            
            correct = np.sum(pred_val - current_label[start_idx:end_idx] <= 0.001)
            print("correct: ", correct)
            #correct = np.sum(pred_val == current_label[start_idx:end_idx])
            total_correct += correct
            total_seen += cur_batch_size
            loss_sum += batch_loss_sum
            
            for i in range(start_idx, end_idx):
                
                l = current_label[i]
                #total_seen_class[l] += 1
                #total_correct_class[l] += (pred_val[i-start_idx] - l <= 0.001)
                fout.write('%d, %d\n' % (pred_val[i-start_idx], l))
                
                if pred_val[i-start_idx] - l > 0.001 and FLAGS.visu: # ERROR CASE, DUMP!
                    
                    out_file = open(out_path, 'a')
                    out_file.write("%s\n" % str(current_data[i].tolist()))
                    out_file.close()
                    
                    '''
                    img_filename = '%d_label_%s_pred_%s.jpg' % (error_cnt, SHAPE_NAMES[l],
                                                           SHAPE_NAMES[pred_val[i-start_idx]])
                    img_filename = os.path.join(DUMP_DIR, img_filename)
                    output_img = pc_util.point_cloud_three_views(np.squeeze(current_data[i, :, :]))
                    scipy.misc.imsave(img_filename, output_img)
                    error_cnt += 1
                    '''
                
    log_string('eval mean loss: %f' % (loss_sum / float(total_seen)))
    log_string('eval accuracy: %f' % (total_correct / float(total_seen)))
    '''
    log_string('eval avg class acc: %f' % (np.mean(np.array(total_correct_class)/np.array(total_seen_class,dtype=np.float))))
    
    class_accuracies = np.array(total_correct_class)/np.array(total_seen_class,dtype=np.float)
    for i, name in enumerate(SHAPE_NAMES):
        log_string('%10s:\t%0.3f' % (name, class_accuracies[i]))
    '''


if __name__=='__main__':
    with tf.Graph().as_default():
        evaluate(num_votes=1)
    LOG_FOUT.close()
