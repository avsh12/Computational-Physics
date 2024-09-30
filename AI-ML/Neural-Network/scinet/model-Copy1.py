import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import Input, Model
from tensorflow.keras.layers import Dense

#dfine a sampling subclass of Layer
#the Sampling layer takes two inputs: mean and log_sigma
#the layer samples random samples from mean and log_sigma assuming that the mean and log_sigma are from from a gaussian distribution
class Sampling(keras.layers.Layer):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
    def call(self, inputs):
        #unpack the inputs
        mean, log_var = inputs
        batch_size = tf.shape(mean)[0]
        #size of the latent layer
        units = tf.shape(mean)[1]
        #generate random samples from a gaussian distribution of mean=0 and standard deviation=1
        epsilon = tf.random.normal(shape=(batch_size, units), mean=0, stddev=1)
        #shift and scale epsilon
        return mean + tf.math.exp(log_var/2)*epsilon

class Evolution(keras.layers.Layer):
    def __init__(self, latent_size=3, **kwargs):
        super().__init__(**kwargs)
        #self.latent_size = latent_size

   def build(self, input_shape):
       self.input_dim = input_shape[-1]

       self.biases = self.add_weights(shape=(self.input_dim,), initializer='random_normal')

   def call(self, inputs):
       return inputs + self.biases

class Encoder(keras.layers.Layer):
    def __init__(self, latent_size=3, layer_size=[70, 70], activation='elu', vae=True, **kwargs):
        super().__init__(**kwargs)
        self.layer_size = layer_size
        self.latent_size = latent_size
        self.vae = vae
        
        self.encoder_layer = []

        for i in range(len(self.layer_size)):
            self.encoder_layer.append(Dense(layer_size[i], activation=activation, name=('encoder_layer'+str(i+1))))

        if self.vae:
            self.mean = Dense(latent_size, name='mean')
            self.log_var = Dense(latent_size, name='log_sigma')
            #generate random samples from gaussian distribution with mean = mean and log(variance) = log_var
            self.z = Sampling(name='sampling_layer')
        else:
            self.z = Dense(latent_size, name='latent_layer')
    def call(self, inputs):
        x = inputs
        for i in range(len(self.layer_size)):
            x = self.encoder_layer[i](x)
        if self.vae:
            mean = self.mean(x)
            log_var = self.log_var(x)
            z = self.z([mean, log_var])
            return mean, log_var, z
        else:
            z = self.z(x)
            return z

class Decoder(keras.layers.Layer):
    def __init__(self, layer_size=[70, 70], output_size=1, activation='elu', **kwargs):
        super().__init__(**kwargs)
        self.layer_size = layer_size
        self.output_size = output_size
        
        self.decoder_layer = []
        for i in range(len(self.layer_size)):
            self.decoder_layer.append(Dense(layer_size[i], activation=activation, name=('encoder_layer'+str(i+1))))

        self.decoder_output = Dense(output_size, name='decoder_layer_out')
        
    def call(self, inputs):
        x = inputs
        for i in range(len(self.layer_size)):
            x = self.decoder_layer[i](x)
            
        x = self.decoder_output(x)
        return x

#create a subclass VAE of Model
#the model takes encoder and decoder models and joins them to build beta-VAE
#the model defines a training method to train the beta-VAE model so constructed
"""Variables in the class object:

   encoder               --   encoder model to produce mean and log_var as output from the input
   decoder               --   decoder model to produce the output
   reconstruction_loss   --   total loss for the output accuracy of the decoder - sum over total outputs and average over samples
   kl_loss               --   total loss for the kl-metric - sum over output units of latent layer and average over samples
   total_loss            --   reconstruction_loss + kl_loss
   beta_rec              --   parameter to control the relative significance of reconstruction loss
   beta_kl               --   parameter to control the relative significance of kl-divergence loss
"""

class ScNN(Model):
    #
    def __init__(self, input_size=50, encoder_layer=[70, 70], latent_size=3, question_size=1, decoder_layer=[70, 70], output_size=1, beta_rec=500., beta_kl=1., target_var=0.01, activation='elu', **kwargs):
        super().__init__(**kwargs)

        self.encoder = Encoder(layer_size=encoder_layer, latent_size=latent_size, activation=activation, vae=True, name='encoder')
        if self.question_size != 0:
            self.concatenate = keras.layers.Concatenate()
        self.decoder = Decoder(layer_size=decoder_layer, output_size=output_size, activation=activation, name='decoder')

        self.question_size = question_size
        self.beta_rec = beta_rec
        self.beta_kl = beta_kl
        self.target_var = target_var

        #define the losses
        #Mean stores the mean value of the argument passed to it and a name for the mean value
        self.total_loss = keras.metrics.Mean(name='total_loss')
        self.total_reconstruction_loss = keras.metrics.Mean(name='reconstruction_loss')
        self.total_kl_loss = keras.metrics.Mean(name='kl_loss')

    #return the losses to be printed during training
    @property
    def metrics(self):
        return [self.total_loss, self.total_reconstruction_loss, self.total_kl_loss]

    #feeds an output of encoder and question to the decoder and returns the outputs of encoder and decoder
    #encoder input and a question is passed as input
    def call(self, inputs):
        if self.question_size != 0:
            observation, question = inputs
        else:
            observation = inputs
            
        mean, log_var, z = self.encoder(observation)

        if self.question != 0:
            decoder_input = self.concatenate([z, question])
            predictions = self.decoder(decoder_input)
        else:
            predictions self.decoder(z)
            
        return mean, log_var, z, predictions

    #training method called when called to fit
    #the method is passed a set of observation, a question for each observation and the corresponding answers
    def train_step(self, data):
        inputs = data['time_series']
        if self.question_size != 0:
            question = data['t_question']
        outputs = data['answer']

        with tf.GradientTape() as tape:
            #feed-forward step
            #calls the call() method
            mean, log_var, z, predictions = self([inputs, question]) if self.question_size != 0 else self(inputs)

            #calculate losses for the current feed-forward step
            #outputs and predictions are one dimensional so only take mean over samples
            reconstruction_loss = self.beta_rec*tf.reduce_mean(tf.reduce_sum((outputs - predictions)**2, axis=1))
            #take sum over dimension and mean over samples
            #kl_loss = 0.5*self.beta_kl*tf.reduce_mean(tf.reduce_sum(tf.square(z_mean) + 2*tf.exp(z_log_sigma) -  2*z_log_sigma, axis=1))
            kl_loss = 0.5*self.beta_kl*tf.reduce_mean(tf.reduce_sum(tf.square(mean)/self.target_var + tf.exp(log_var)/self.target_var - log_var + tf.math.log(self.target_var), axis=1)-3)
            loss = reconstruction_loss + kl_loss

            #calculates gradients w.r.t. the trainable weights of the model
            gradients = tape.gradient(loss, self.trainable_weights)
            #updates the trainable weights using the gradients calculated
            #method to update the weights is passed in call to compile function of the model as optimizer parameter
            self.optimizer.apply_gradients(zip(gradients, self.trainable_weights))

            #updates the current losses
            self.total_reconstruction_loss.update_state(reconstruction_loss)
            self.total_kl_loss.update_state(kl_loss)
            self.total_loss.update_state(loss)

            #return the losses name and the corresponding values
            return {m.name: m.result() for m in self.metrics}