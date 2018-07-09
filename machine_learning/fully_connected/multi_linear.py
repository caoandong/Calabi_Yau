import tensorflow as tf
import numpy as np
import math
import sys
import os

# Input: parametrized heights
# Output: volume

def placeholder_inputs(batch_size, num_param):
    param_pl = tf.placeholder(name="input_param", dtype=tf.float32, shape=(batch_size, num_param))
    vol_pl = tf.placeholder(name="input_vol", dtype=tf.float32, shape=(batch_size))
    
    return param_pl, vol_pl

def load_input(train_path):
    train = open(train_path, 'r')
    train_x = []
    train_y = []
    for line in train:
        data = eval(line)
        if data[1] > 0:
            train_x.append(data[0])
            inv_vol = 1/data[1]
            train_y.append(inv_vol)
    print ('Number of data: ', len(train_x))
    
    return train_x, train_y

def expand_2(input_param_1, batch_size, num_param):
    for row in range(batch_size):
        param_tmp = input_param_1[row,:]
        param_tmp = tf.reshape(param_tmp, [1, -1])
        param_mat = tf.matmul(tf.transpose(param_tmp),param_tmp)

        for i in range(num_param):
            param_tmp2 = tf.reshape(param_mat[i,i:], [1,-1])
            param_tmp = tf.concat([param_tmp,param_tmp2], 1)
        param_tmp = tf.reshape(param_tmp, [-1])

        try:
            input_parameter = tf.concat([input_parameter,param_tmp[None, :]],0)
            input_parameter = tf.reshape(input_parameter, [row+1, -1])
        except:
            input_parameter = param_tmp
            input_parameter = tf.reshape(input_parameter, [row+1, -1])
    
    return input_parameter

def fl_linear(param, batch_size, num_param):
    with tf.variable_scope("fl_linear", reuse=True):
        input_param = tf.placeholder(name="input_param", dtype=tf.float32, shape=(batch_size, num_param))
        input_vol = tf.placeholder(name="input_vol", dtype=tf.float32, shape=(batch_size))
        # Linear activation
        net = tf.fully_connected(1, activation_fn=None, scope='fc1')
    
    return net

def mat_mul(param, batch_size, num_param):
    with tf.variable_scope("mat_mul", reuse=True):
        input_param_1 = tf.placeholder(name="input_param_1", dtype=tf.float32, shape=(batch_size, num_param))
        input_volume = tf.placeholder(name="input_volume", dtype=tf.float32, shape=(batch_size))
        
        input_parameter = expand_2(input_param_1, batch_size, num_param)
        num_param = int(num_param + num_param*(num_param+1)/2)

        # Layer 1:
        # Input: (batch_size, num_param*2)
        # W1: (num_param*2, 1), b: (batch_size, 1)
        # Output: (batch_size, 1)
        W1 = tf.Variable(tf.random_normal([num_param, 1], dtype=tf.float32), name="W")
        b1 = tf.Variable(tf.zeros([batch_size, 1], dtype=tf.float32), name="b")
        y = tf.add(tf.matmul(input_parameter, W1), b1)
        
        cost_op = tf.reduce_mean(tf.pow(y-input_volume, 2))
    
    return y, cost_op, input_param_1, input_volume, input_parameter

def mat_mul_2(param, batch_size, num_param):
    with tf.variable_scope("mat_mul_2", reuse=True):
        input_parameter = tf.placeholder(name="input_parameter", dtype=tf.float32, shape=(batch_size, num_param))
        input_volume = tf.placeholder(name="input_volume", dtype=tf.float32, shape=(batch_size))
        
        # Layer 1:
        # Input: (batch_size, num_param)
        # W1: (num_param, 10), b: (batch_size, 10)
        # Output: (batch_size, 10)
        W1 = tf.Variable(tf.random_normal([num_param, 10], dtype=tf.float32), name="W1")
        b1 = tf.Variable(tf.zeros([batch_size, 10], dtype=tf.float32), name="b1")
        y = tf.add(tf.matmul(input_parameter, W1), b1)
        
        # Layer 2:
        # Input: (batch_size, 10)
        # W1: (10, 1), b: (batch_size, 1)
        # Output: (batch_size, 1)
        W2 = tf.Variable(tf.random_normal([10, 1], dtype=tf.float32), name="W2")
        b2 = tf.Variable(tf.zeros([batch_size, 1], dtype=tf.float32), name="b2")
        y = tf.add(tf.matmul(y, W2), b2)
        
        cost_op = tf.reduce_mean(tf.pow(y-input_volume, 2))

    return y, cost_op, input_parameter, input_volume

def train(data, vol):
    num_data = len(data)
    split_idx = int(num_data*0.8)
    
    input_data = data[0:split_idx]
    input_vol = vol[0:split_idx]
    
    test_data = data[split_idx:]
    test_vol = vol[split_idx:]
    
    num_train = split_idx
    num_test = num_data - split_idx
    num_param = len(input_data[0])
    num_test = len(test_data)
    print('number of train: ', num_train)
    print('number of test: ', num_test)
    print('number of parameters: ', num_param)
    print('number of test: ', num_test)
    
    num_batches = int(num_train/batch_size)
    test_batches = int(num_test/batch_size)
    input_param = np.array(input_data)
    input_vol = np.array(input_vol)
    test_param = np.array(test_data)
    test_vol = np.array(test_vol)
    
    #y = fl_linear(input_param)
    y, cost_op, input_parameter, input_volume, input_param_concat = mat_mul(input_param, batch_size, num_param)
    
    learning_rate = tf.Variable(0.5, trainable=False)
    optimizer = tf.train.AdamOptimizer(learning_rate)
    train_op = optimizer.minimize(cost_op)
    
    sess = tf.Session() # Create TensorFlow session
    saver = tf.train.Saver() # Save Tensorflow graph
    print ("Beginning Training")
    with sess.as_default():
        init = tf.global_variables_initializer()
        sess.run(init)
        sess.run(tf.assign(learning_rate, alpha))
        for epoch in range(1, max_epochs):
            print('epoch: ', epoch)
            train_one_epoch(sess, cost_op, train_op, num_batches, input_parameter, input_param, input_volume, input_vol, input_param_concat, y)
            eval_one_epoch(sess, cost_op, test_batches, input_parameter, test_param, input_volume, test_vol)
            
            
            
            if epoch % 10 == 0:
                path = "log/model.ckpt"
                save_path = saver.save(sess, path)
                print("Saved.")
            
        print('Done.')
        
        collection = tf.get_collection(tf.GraphKeys.GLOBAL_VARIABLES, scope='mat_mul')
        for x in collection:
            if x.name == 'mat_mul/W:0' or x.name == 'mat_mul/b:0':
                print("name: ", x.name)
                print("variables: ", sess.run(x))
            

def train_one_epoch(sess, cost_op, train_op, num_batches, input_parameter, input_param, input_volume, input_vol, input_param_concat, y):
     for batch_idx in range(num_batches):
        start_idx = batch_idx * batch_size
        end_idx = (batch_idx+1) * batch_size
        
        cost,_ = sess.run([cost_op, train_op], feed_dict={input_parameter: input_param[start_idx:end_idx], input_volume: input_vol[start_idx:end_idx]})
        
        #concat, pred = sess.run([input_param_concat, y], feed_dict={input_parameter: input_param[start_idx:end_idx], input_volume: input_vol[start_idx:end_idx]})
        #print("input volume: ", input_vol[start_idx:end_idx])
        #print("predicted value: ", pred)
        '''
        try:
            cost,_ = sess.run([cost_op, train_op], feed_dict={input_parameter: input_param[start_idx:end_idx], input_volume: input_vol[start_idx:end_idx]})
            
            
            #print("input parameter: ", input_param[start_idx:end_idx])
            #concat, pred = sess.run([input_param_concat, y], feed_dict={input_parameter: input_param[start_idx:end_idx], input_volume: input_vol[start_idx:end_idx]})
            #print("parameter concatinated: ", concat)
            print("input volume: ", input_vol[start_idx:end_idx])
            print("predicted value: ", pred)
            
        except:
            print("input error: ")
            print("input parameter: ", input_param[start_idx:end_idx])
            print("input volume: ", input_volume[start_idx:end_idx])
        '''
        print("loss: ", cost)
        
def eval_one_epoch(sess, cost_op, test_batches, input_parameter, test_param, input_volume, test_vol):
    for batch_idx in range(test_batches):
        start_idx = batch_idx * batch_size
        end_idx = (batch_idx+1) * batch_size
        cost = sess.run(cost_op, feed_dict={input_parameter: test_param[start_idx:end_idx], input_volume: test_vol[start_idx:end_idx]})
        print("test loss: ", cost)

# Parameters:
last_cost = 0
alpha = 0.4
batch_size = 10
max_epochs = 10000
tolerance = 1e-3

    
train_path = '/home/carnd/CYML/output/train/cylinder/tri_1_to_50_2.txt'
train_data, train_vol = load_input(train_path)
print ('sample data: ', train_data[0])
train(train_data, train_vol)