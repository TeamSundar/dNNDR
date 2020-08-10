from keras.models import Model
from keras.regularizers import l2
from keras.constraints import max_norm
from keras.layers import Input, Dense, Dropout, Flatten, Activation, Concatenate, Layer
from keras.layers import Conv1D, Add, MaxPooling1D, BatchNormalization
from keras.layers import Embedding, Bidirectional, GlobalMaxPooling1D, LSTM, CuDNNLSTM
import keras.backend as K

from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession

# Attention class
class attention(Layer):
    def __init__(self,**kwargs):
        super(attention,self).__init__(**kwargs)

    def build(self,input_shape):
        self.W=self.add_weight(name="att_weight",shape=(input_shape[-1],1),initializer="normal")
        self.b=self.add_weight(name="att_bias",shape=(input_shape[1],1),initializer="zeros")        
        super(attention, self).build(input_shape)

    def call(self,x):
        et=K.squeeze(K.tanh(K.dot(x,self.W)+self.b),axis=-1)
        at=K.softmax(et)
        at=K.expand_dims(at,axis=-1)
        output=x*at
        return K.sum(output,axis=1)

    def compute_output_shape(self,input_shape):
        return (input_shape[0],input_shape[-1])

    def get_config(self):
        return super(attention,self).get_config()

class ARCHITECTURE():
    def __init__(self):
        pass
    
    def protSeq(self, input_shape):
        input_target = Input(shape=(input_shape,))
        emb_target = Embedding(21, 128, input_length=1000)(input_target) 
        conv_target_1 = Conv1D(filters=32, kernel_size=3, padding='same', activation='relu')(emb_target)
        pool_target_1 = MaxPooling1D(pool_size=2)(conv_target_1)
        att_in_target = Bidirectional(CuDNNLSTM(32, kernel_regularizer=l2(0.01), return_sequences=True, recurrent_regularizer=l2(0.01), bias_regularizer=l2(0.01)))(pool_target_1)
        att_out_target = attention()(att_in_target) 
        return input_target, att_out_target

    def drugSeq(self, input_shape):
        input_drug = Input(shape=(input_shape,))
        emb_drug = Embedding(44, 128, input_length=1000)(input_drug) 
        conv_drug_1 = Conv1D(filters=32, kernel_size=3, padding='same', activation='relu')(emb_drug)
        pool_drug_1 = MaxPooling1D(pool_size=2)(conv_drug_1)
        att_in_drug = Bidirectional(CuDNNLSTM(32, kernel_regularizer=l2(0.02), return_sequences=True, recurrent_regularizer=l2(0.02), bias_regularizer=l2(0.02)))(pool_drug_1)
        att_out_drug = attention()(att_in_drug)
        return input_drug, att_out_drug

    def drugDes(self, input_shape):
        input_drug_des = Input(shape=(input_shape,))
        dense_drug_des_1 = Dense(512, activation="relu", kernel_initializer='glorot_normal')(input_drug_des)
        dense_drug_des_2 = Dense(512, activation="relu", kernel_initializer='glorot_normal')(dense_drug_des_1)
        return input_drug_des, dense_drug_des_2

    def protDes(self, input_shape):
        pass
    
    def drugXAE(self, input_shape):
        pass
    

