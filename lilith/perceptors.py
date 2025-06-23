import files
import os
import sys



class MLPerceptron:
	def __init__(self, epochs, input_size, learnrate):
		self.epochs = epochs
		self.input_size = input_size
		self.learnrate = learnrate
		self.weights = None #i think numpy would be good, time to learn venv
		self.bias = None #again, numpy for initial weights?
		self.misses = []

	def predict(self, ins):


	def activator():