
import sys

import warnings
warnings.filterwarnings("ignore")

import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation
from keras.optimizers import SGD, RMSprop, Adam, Nadam
from keras.callbacks import Callback, EarlyStopping, ModelCheckpoint
from keras import backend, regularizers
from keras.layers.advanced_activations import LeakyReLU
from keras.utils import to_categorical, normalize

import numpy as np
import pandas as pd
import time, os

from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import LabelEncoder

%matplotlib inline
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from IPython import display

use_gpu = True
if use_gpu:
    import tensorflow as tf
    config = tf.ConfigProto( device_count = {'GPU': 1 , 'CPU': 8} ) 
    sess = tf.Session(config=config) 
    backend.set_session(sess)
    
%load_ext autoreload
%autoreload 2