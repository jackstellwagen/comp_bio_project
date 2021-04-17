import tensorflow as tf

from tensorflow.keras import datasets, layers, models
import numpy as np
from seq_dataset import SeqDataSet
from toolkit import ReadDataSet



stage_depth = 1#10 # the number of BasicBlocks in a Stage (N)
kernel_width = 5 # Height (k)
amino_dim = 30 #the number of features, NOT the number of amino acids
std_in_channel = 256 # (c)
std_out_channel = std_in_channel * 2

std_filter_shape = [kernel_width, 1, std_in_channel, std_out_channel]
attn_filter_shape = [kernel_width, 1, std_in_channel, std_in_channel]

filter_shape = [3, amino_dim, 1, std_out_channel]




def add_GLU(input):
    #split the input in to 2 tensors
    half1, half2 = layers.Lambda( lambda x: tf.split(x,num_or_size_splits=2,axis=3))(input)

    #multiply the two halves together and return
    # This is how GLU is defined line 49 tf_lib
    return layers.Multiply()([half1,(layers.Activation("sigmoid")(half1))])




def add_NormBlock(input):

    #might need to specify axis TODO
    normalized_input = layers.LayerNormalization()(input)

    GLU = add_GLU(normalized_input)

    return GLU

def add_ResidualBlock(input):
    padding = std_filter_shape[0] // 2

    GLU = add_NormBlock(input)

    padded_GLU = layers.ZeroPadding2D(padding = ((padding,padding), (0,0)) ) (GLU)

    output = layers.Conv2D(std_out_channel, (kernel_width,1),use_bias = False) (padded_GLU)

    #add the output and the input
    #thats what makes this "residual"
    #return output
    #sizes dont add up here for some reason, this should return below but I need
    # to find the bug
    return layers.Add()([input,output])



def add_PlainBlock(input):
    padding = std_filter_shape[0] // 2

    GLU = add_NormBlock(input)

    padded_GLU = layers.ZeroPadding2D(padding = ((padding,padding), (0,0)) ) (GLU)

    output = layers.Conv2D(std_out_channel, (kernel_width,1),use_bias = False) (padded_GLU)

    return output





def get_model():

    # In figure 2 this is the first image to the second image

    input = layers.Input(shape = (None, amino_dim, 1))
    padded_input = layers.ZeroPadding2D(padding = ((1,1), (0,0)) ) (input)
    conv0 = layers.Conv2D(std_out_channel, (3,amino_dim),use_bias = False ) (padded_input)

    #Stage 1, 3rd image
    # stage_depth ResidualBlocks are added followed by 1 plain block
    stage1 = conv0
    for i in range(stage_depth):
        stage1 = add_ResidualBlock(stage1)
    stage1 = add_PlainBlock(stage1)

    #Stage 2
    # stage_depth ResidualBlocks are added followed by 1  NormBlock
    stage2 = stage1
    for i in range(stage_depth):
        stage2 = add_ResidualBlock(stage2)
    stage2 = add_NormBlock(stage2)


    #Final stage, called "proj" in the code

    #this is the part that goes from stage 2 to the final image in fig 2
    conv1 = layers.Conv2D(std_in_channel, (1,1), activation='relu' ) (stage2)

    # Reference code has a dropout layer with a trainable parameter here
    # dropout = layers.Dropout()

    #this is the last convolution in fig 2
    conv2 = layers.Conv2D(2, (1,1), activation = 'relu') (conv1)

    #they apply softmax at the end of the network
    #this could probably be done instead of relu above
    output = layers.Activation("softmax")(conv2)
    model = tf.keras.Model(input, output)


    return model


def loss_fn(prediction, label):
    length = label.shape[1]
    print(label.shape, prediction.shape)
    # prediction = tf.argmax(prediction,1)
    # label = tf.argmax(label,1)

    prediction_slice = prediction[:,:length,:]
    print(prediction_slice.shape)
    return tf.keras.losses.MSE(prediction_slice, label)


def get_raw_data(dataset_name):
    dataset_dir = "DataSet/" + dataset_name + "/"

    dataset = ReadDataSet(dataset_dir, 'data.feat', 'data.lab', SeqDataSet, dataset_name)
    dataset.SetSignal()
    batch_x, batch_y, lens_x, lens_y = dataset.NextRestrictedPaddingBatch(10000)
    x = np.array( [np.array(x) for x in batch_x])

    x_data = np.reshape(np.array(batch_x), newshape = (x.shape[0], x.shape[1]//30, 30))
    #print(x_data)

    y = np.array(batch_y)
    y_data = np.reshape(y, newshape=(y.shape[0], y.shape[1], 1, y.shape[2]))
    return x_data, y_data


def load_data():



    SITA_x, SITA_y = get_raw_data("SITA")
    SITA_EX1_x, SITA_EX1_y = get_raw_data("SITA_EX1")
    SITA_EX2_x,SITA_EX2_y = get_raw_data("SITA_EX2")
    SITA_EX3_x, SITA_EX3_y = get_raw_data("SITA_EX3")

    print(SITA_x.shape)
    print(SITA_EX1_x.shape)
    print(SITA_EX2_x.shape)
    print(SITA_EX3_x.shape)

    train_input = np.concatenate((SITA_EX3_x,SITA_EX1_x))
    train_input = np.concatenate((train_input, SITA_EX2_x))

    train_labels = np.concatenate((SITA_EX3_y,SITA_EX1_y))
    train_labels = np.concatenate((train_labels, SITA_EX2_y))




    return train_input, train_labels, SITA_x, SITA_y

inputs, labels, test_inputs, test_labels = load_data()



model = get_model()
model.compile(optimizer = "adam", loss = loss_fn)
model.summary()

model.fit(
    x=inputs,
    y=labels,
    epochs=2,
    validation_data=(test_inputs, test_labels),
    batch_size=12,
)
