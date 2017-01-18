# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 16:49:34 2016

@author: Ashiqul
"""
from __future__ import division
import pandas as pd
import numpy as np
from sklearn import datasets 
from sklearn.cross_validation import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Perceptron
from sklearn.metrics import accuracy_score
from sklearn.feature_extraction.text import CountVectorizer
from Bio import SeqIO
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_score
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier 
from feature_extraction import extract_motifs_pos 
from sklearn.model_selection import ShuffleSplit


def classifier(maxWordLength, maxDistLength, wordMotifMaxBuffer, distMotifMaxBuffer):
	corpus1=[]#make a blank list of corpus for the training set
	tag=[]#make the list of outcomes
	true = open("ch.fasta")
	false = open("nonch.fasta")
	for line in SeqIO.parse(true, "fasta"):
	    line = line.seq.tostring().lower().replace("x","")
	    line = line.replace('-', "")
	   # print line
	    tag.append("1")
	    fullstring = extract_motifs_pos(line, 1, maxWordLength, 1, maxDistLength, wordMotifMaxBuffer, distMotifMaxBuffer)
	    #fullstring = fullstring+ " "+ pos_prot_1st_word(line)
	    corpus1.append(fullstring) #apperd string from each protein to corpus
	true.close()
	for line in SeqIO.parse(false, "fasta"):
	    line = line.seq.tostring().lower().replace("x","")
	    line = line.replace('-', "")
	    #print line
	    tag.append("0")
	    fullstring = extract_motifs_pos(line, 1, maxWordLength, 1, maxDistLength, wordMotifMaxBuffer, distMotifMaxBuffer)
	    #fullstring = fullstring+ " "+ pos_prot_1st_word(line)
	    corpus1.append(fullstring) #apperd string from each protein to corpus 
	false.close()
	
	corpus = np.array(corpus1) #convert corpus into numpy array
	tag = np.array(tag)  # convert tag into numpy array   
	#print corpus # print for debugging 
	#print tag # print for debugging
	
	count = CountVectorizer(max_features=15000000, vocabulary = None, max_df=0.3, min_df = 3, stop_words=[1,2])#giving the CountVectorizer function a short name
	#get the vocabulary of train set to use for the test set
	
	
	
	bag = count.fit_transform(corpus) #transform the corpus(bag of words into sparse martrix)
	#print (count.vocabulary_) #count the occurence of number of words
	##get the vocabulary of train set to use for the test set. Next time put the "voc" in
	#the vocabulary parameter of count
	voc = count.get_feature_names() 
	#print len(voc)
	bag= bag.toarray() #convert the sparsematrix into an array
	np.place(bag, bag>0, [1])
	#print bag
	
	forest = RandomForestClassifier(n_estimators = 1000,
	                                random_state = 1,
	                                n_jobs =1)
	forest.fit(bag[:, 0:-1], tag)
	importances = forest.feature_importances_
	std = np.std([tree.feature_importances_ for tree in forest.estimators_],
	             axis=0)
	indices = np.argsort(importances)[::-1]
	
	# Print the feature ranking
	#print("Feature ranking:")
	important = list()
	for f in range(0,500):
	    important.append(indices[f])
	bag=bag[:,important]
	bag=pd.DataFrame(bag)
	bag['tag']=tag
	voc = np.array(voc)[important]
	x = bag.iloc[:, 0:-1]
	y = bag.iloc[:,-1]
	#parameterize the Logistic Regression algorithm
	cv = ShuffleSplit(n_splits=10, test_size=0.2)
	clf= MLPClassifier(solver='lbfgs', alpha=1e-5,hidden_layer_sizes=(15), random_state=1)
	return np.mean(cross_val_score(clf, x, y, cv =cv)) 

def main():
	score = []
	for maxWordLength in range(2, 40):
		score.append(classifier(maxWordLength, 15, 1, 0))
	print("maxWordLength 2 - 40:") 
	print(score)
	score = []
	for maxDistLength in range(2,40):
		score.append(classifier(15, maxDistLength, 1, 0))
	print("maxDistLength 2-40")
	print(score)
	score = []
	for wordMotifMaxBuffer in range(0, 10):
		score.append(classifier(15, 15, wordMotifMaxBuffer, 0))
	print("wordMotifMaxBuffer 0-10")
	print(score)
	score = []
	for distMotifMaxBuffer in range(0, 10):
		score.append(classifier(15, 15, 1, distMotifMaxBuffer))
	print("distMotifMaxBuffer 0-10")	
	print(score)

if __name__ == "__main__": 
	main()	
