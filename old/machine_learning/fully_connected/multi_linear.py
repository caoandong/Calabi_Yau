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

def load_input_extend(train_path):
    train = open(train_path, 'r')
    train_x = []
    train_y = []
    for line in train:
        data = eval(line)
        if data[1] > 0:
            p1 = data[0][0]
            p2 = data[0][1]
            p3 = data[0][2]
            p = p1**2
            data[0].append(p)
            p = p2**2
            data[0].append(p)
            p = p3**2
            data[0].append(p)
            p = p1*p2
            data[0].append(p)
            p = p1*p3
            data[0].append(p)
            p = p2*p3
            data[0].append(p)
            p = p1**3
            data[0].append(p)
            p = p2**3
            data[0].append(p)
            p = p3**3
            data[0].append(p)
            p = (p1**2)*p2
            data[0].append(p)
            p = (p1**2)*p3
            data[0].append(p)
            p = p1*(p2**2)
            data[0].append(p)
            p = p1*(p3**2)
            data[0].append(p)
            p = (p2**2)*p3
            data[0].append(p)
            p = (p3**2)*p2
            data[0].append(p)
            p = p1*p2*p3 
            data[0].append(p)
            train_x.append(data[0])
            inv_vol = float(data[1])**(-1)
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
    
    num_param = int(num_param + num_param*(num_param+1)/2)
    
    return input_parameter, num_param

def expand_3(input_param_1, batch_size, num_param):
    # TODO: make this more general
    num_param = 3
    input_param = 10
    for row in range(batch_size):
        param_tmp = input_param_1[row,:]
        param_tmp = tf.reshape(param_tmp, [1, -1])
        p1 = param_tmp[0,0]
        p2 = param_tmp[0,1]
        p3 = param_tmp[0,2]
        
        p = tf.pow(p1,2)
        p = tf.reshape(p, [-1])
        param_tmp = tf.concat([param_tmp, p[None,:]], 1)
        p = tf.pow(p2,2)
        p = tf.reshape(p, [-1])
        param_tmp = tf.concat([param_tmp, p[None,:]], 1)
        p = tf.pow(p3,2)
        p = tf.reshape(p, [-1])
        param_tmp = tf.concat([param_tmp, p[None,:]], 1)
        
        p = tf.multiply(p1, p2)
        p = tf.reshape(p, [-1])
        param_tmp = tf.concat([param_tmp, p[None,:]], 1)
        p = tf.multiply(p1, p3)
        p = tf.reshape(p, [-1])
        param_tmp = tf.concat([param_tmp, p[None,:]], 1)
        p = tf.multiply(p2, p3)
        p = tf.reshape(p, [-1])
        param_tmp = tf.concat([param_tmp, p[None,:]], 1)
        
        p = tf.pow(p1,3)
        p = tf.reshape(p, [-1])
        param_tmp = tf.concat([param_tmp, p[None,:]], 1)
        p = tf.pow(p2,3)
        p = tf.reshape(p, [-1])
        param_tmp = tf.concat([param_tmp, p[None,:]], 1)
        p = tf.pow(p3,3)
        p = tf.reshape(p, [-1])
        param_tmp = tf.concat([param_tmp, p[None,:]], 1)
        
        p = tf.multiply(tf.pow(p1,2), p2)
        p = tf.reshape(p, [-1])
        param_tmp = tf.concat([param_tmp, p[None,:]], 1)
        p = tf.multiply(tf.pow(p1,2), p3)
        p = tf.reshape(p, [-1])
        param_tmp = tf.concat([param_tmp, p[None,:]], 1)
        p = tf.multiply(p1, tf.pow(p2,2))
        p = tf.reshape(p, [-1])
        param_tmp = tf.concat([param_tmp, p[None,:]], 1)
        p = tf.multiply(p1, tf.pow(p3,2))
        p = tf.reshape(p, [-1])
        param_tmp = tf.concat([param_tmp, p[None,:]], 1)
        
        p = tf.multiply(tf.pow(p2,2), p3)
        p = tf.reshape(p, [-1])
        param_tmp = tf.concat([param_tmp, p[None,:]], 1)
        p = tf.multiply(tf.pow(p3,2), p2)
        p = tf.reshape(p, [-1])
        param_tmp = tf.concat([param_tmp, p[None,:]], 1)
        
        p = tf.multiply(tf.multiply(p1, p2),p3)
        p = tf.reshape(p, [-1])
        param_tmp = tf.concat([param_tmp, p[None,:]], 1)
        
        try:
            param_tmp = tf.reshape(param_tmp, [-1])
            input_param = tf.concat([input_param, param_tmp[None,:]], 0)
        except:
            input_param = tf.reshape(param_tmp, [1, -1])
    
    num_param = num_param + 6+ 10
    return input_param, num_param

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
        
        #input_parameter, num_param = expand_2(input_param_1, batch_size, num_param)
        input_parameter, num_param = expand_3(input_param_1, batch_size, num_param)
        
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
        
        #cost_op = tf.reduce_mean(tf.multiply(tf.pow(tf.pow(y,-1)-tf.pow(input_volume,-1), 2), y))
        cost_op = tf.reduce_mean(tf.pow(y-input_volume, 2))

    return y, cost_op, input_parameter, input_volume

def mat_mul_3(param, batch_size, num_param):
    # TODO: make this more general
    with tf.variable_scope("mat_mul_3", reuse=True):
        input_parameter = tf.placeholder(name="input_parameter", dtype=tf.float32, shape=(batch_size, num_param))
        input_volume = tf.placeholder(name="input_volume", dtype=tf.float32, shape=(batch_size))
        
        # Layer 1:
        # Input: (batch_size, num_param)
        # W1: (num_param, 10), b: (batch_size, 10)
        # Output: (batch_size, 10)
        W1 = tf.Variable(tf.random_normal([num_param, 10], dtype=tf.float32), name="W1")
        y = tf.matmul(input_parameter, W1)
        
        # Layer 2:
        # Input: (batch_size, 10)
        # W1: (10, 5), b: (batch_size, 5)
        # Output: (batch_size, 5)
        W2 = tf.Variable(tf.random_normal([10, 5], dtype=tf.float32), name="W2")
        y = tf.matmul(y, W2)
        
        # Layer 3:
        # Input: (batch_size, 5)
        # W1: (5, 1), b: (batch_size, 1)
        # Output: (batch_size, 1)
        W3 = tf.Variable(tf.random_normal([5, 1], dtype=tf.float32), name="W3")
        b3 = tf.Variable(tf.zeros([batch_size, 1], dtype=tf.float32), name="b3")
        y = tf.add(tf.matmul(y, W3), b3)
        
        #cost_op = tf.reduce_mean(tf.pow(tf.pow(y,-1)-tf.multiply(tf.pow(input_volume,-1), y), 2))
        cost_op = tf.reduce_mean(tf.pow(y-input_volume, 2))
        #cost_op = tf.reduce_mean(tf.multiply(tf.pow(y-input_volume, 2), tf.pow(y, -3/4)))

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
    print('batch_size: ', batch_size)
    print('number of batches: ', num_batches)
    test_batches = int(num_test/batch_size)
    input_param = np.array(input_data)
    input_vol = np.array(input_vol)
    test_param = np.array(test_data)
    test_vol = np.array(test_vol)
    
    #y = fl_linear(input_param)
    y, cost_op, input_parameter, input_volume = mat_mul_3(input_param, batch_size, num_param)
    
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
            train_one_epoch(sess, cost_op, train_op, num_batches, input_parameter, input_param, input_volume, input_vol, y)
            eval_one_epoch(sess, cost_op, test_batches, input_parameter, test_param, input_volume, test_vol)
            
            if epoch % 10 == 0:
                path = "log/model_3.ckpt"
                save_path = saver.save(sess, path)
                print("Saved.")
            
        print('Done.')
        
        collection = tf.get_collection(tf.GraphKeys.GLOBAL_VARIABLES, scope='mat_mul_3')
        for x in collection:
            if x.name == 'mat_mul_3/W1:0' or x.name == 'mat_mul_3/W2:0'or x.name == 'mat_mul_3/W3:0' or x.name == 'mat_mul_3/b3:0':
                print("name: ", x.name)
                print("variables: ", sess.run(x))
                output.write("name: %s\n" % x.name)
                output.write("variable: %s\n" % str(sess.run(x)))
            
def train_one_epoch(sess, cost_op, train_op, num_batches, input_parameter, input_param, input_volume, input_vol, y):
    #print('batches: ', num_batches)
    for batch_idx in range(num_batches):
        start_idx = batch_idx * batch_size
        end_idx = (batch_idx+1) * batch_size

        batch_cost,_ = sess.run([cost_op, train_op], feed_dict={input_parameter: input_param[start_idx:end_idx], input_volume: input_vol[start_idx:end_idx]})
        print("loss: ", batch_cost)
        output.write("loss: %f\n" % batch_cost)

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
        
        
def eval_one_epoch(sess, cost_op, test_batches, input_parameter, test_param, input_volume, test_vol):
    for batch_idx in range(test_batches):
        start_idx = batch_idx * batch_size
        end_idx = (batch_idx+1) * batch_size
        batch_cost = sess.run(cost_op, feed_dict={input_parameter: test_param[start_idx:end_idx], input_volume: test_vol[start_idx:end_idx]})
        print("test loss: ", batch_cost)
        output.write("test loss: %f\n" % batch_cost)

# Parameters:
last_cost = 0
alpha = 0.4
batch_size = 10
max_epochs = 10000
tolerance = 1e-3

train_path = '/home/carnd/CYML/output/train/cylinder/lift_1_to_50.txt'
out_path = '/home/carnd/CYML/output/train/cylinder/loss/lift_1_to_50.txt'
output = open(out_path, 'w')
train_data, train_vol = load_input_extend(train_path)
print ('sample data: ', train_data[0])
train(train_data, train_vol)