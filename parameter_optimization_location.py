from __future__ import division
import numpy as np
from sklearn.feature_extraction.text import CountVectorizer
from Bio import SeqIO
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import ShuffleSplit
from sklearn.linear_model import LogisticRegression
from future_module import extract_motifs_pos

fh = open("optimized.txt", "w")
def parameters(f1, f2,  x1, x2, x3, x4, x5, x6, x7, x8):

    corpus1=[]#make a blank list of corpus for the training set
    tag=[]#make the list of outcomes
    true = open("subcellular_localization_without multilabels.fasta")
    count = 1
    for lines in SeqIO.parse(true, "fasta"):
        protid= lines.id
        line = lines.seq.tostring().lower().replace("x","")
        line = line.replace('-', "")
        line = line.upper()
        if "Cytoplasmic" in protid:
            tag.append(0)
        elif "Extracellular" in protid:
            tag.append(1)
        elif "Inner" in protid:
            tag.append(2)
        elif "Outer" in protid:
            tag.append(3)
        elif "Periplasmic" in protid:
            tag.append(4)
        
        
        #print count
        # print line
        fullstring = extract_motifs_pos(line,x1, x2, x3, x4, x5, x6, x7, x8)
        #fullstring = fullstring+ " "+ pos_prot_1st_word(line)
        corpus1.append(fullstring) #apperd string from each protein to corpus
    true.close()
    
    corpus = np.array(corpus1) #convert corpus into numpy array
    tag = np.array(tag)  # convert tag into numpy array   
    #print corpus # print for debugging 
    #print tag # print for debugging
    del fullstring 
    del corpus1
    try:
        count = CountVectorizer(max_features=15000000, vocabulary = None, max_df=0.20, min_df = 7, stop_words=[1,2])#giving the CountVectorizer function a short name
        #get the vocabulary of train set to use for the test set
        
        
        
        bag = count.fit_transform(corpus) #transform the corpus(bag of words into sparse martrix)
        #print (count.vocabulary_) #count the occurence of number of words
        ##get the vocabulary of train set to use for the test set. Next time put the "voc" in
        #the vocabulary parameter of count
        voc = count.get_feature_names() 
        #print len(voc)
        if len(voc)>150000:
            return 0.0
        bag= bag.toarray() #convert the sparsematrix into an array
        np.place(bag, bag>0, [1])
    except:
        return 0.0
    cv = ShuffleSplit(n_splits=5, test_size=0.2, random_state =1)
    clf=  LogisticRegression(C= 1.0, random_state=1)
    return np.mean(cross_val_score(clf, bag, tag, cv =cv)) 



def seeding1(f1, f2, x_seed = 0, iteration = 0, x5=1, x6=1, x7=5, x8=1, step = 3, first_time = True):
   
    for i in range(0,26,step):
            
        #print "seed: ", i
        
        if iteration%2!=0:
            x= i
            x2 = x_seed
            
            #print "x4 is parameter"
            
        else:
            x = i
            x4 = x_seed
            #print "x2 is parameter"
            
        p = 0.0     
        
        if first_time:
            p_new= -1
            p_best = -1
        else:
            p_new = p_best
        first_time = False
        
        used = False
        
        while True:
            if iteration%2!=0:
                x4= x
            else:
                x2 = x
            
            #print "count:  ", count
            #print "value of x4: ", x4
            p = parameters(f1, f2, 0, x2, 0, x4, x5, x6, x7, x8)
            #print "accuracy for seed %d is: %f " % (i, p)
            #print "x4 value is" , x4
            #print "value of p", p
            #print "value of p_new", p_new
            
            if p_new < p:
                p_new = p
                #print "new value of p_new is:  ", p_new
                x+=1
                #print "yes"
                used = True
                #print used
            else:
                #print "want to see if use is true or false  :", used
                #print "no"
              
                p_best = p_new
                if used:
                    x_best = x-1
                    #print 'best accuracy for seed %d is %f' %(i, p_best)
                    #print "best parameter so far: %d" %x_best
                break 
    #print "best accuracy for seed is %f" % p_best
    #print "best parameter so far: %d" %x_best
    return x_best, p_best


def seeding2(f1, f2, x_seed = 0, iteration = 0, x2=1, x4=1, x7=5, x8=1, step = 1, first_time = True):
   
    for i in range(0,11,step):
            
        #print "seed: ", i
        
        if iteration%2!=0:
            x= i
            x5 = x_seed
            
            #print "x4 is parameter"
            
        else:
            x = i
            x6= x_seed
            #print "x2 is parameter"
            
        p = 0.0     
        
        if first_time:
            p_new= -1
            p_best = -1
        else:
            p_new = p_best
        first_time = False
        
        used = False
        
        while True:
            if iteration%2!=0:
                x6= x
            else:
                x5 = x
            
            #print "count:  ", count
            #print "value of x4: ", x4
            p = parameters(f1, f2, 0, x2, 0, x4, x5, x6, x7, x8)
            #print "accuracy for seed %d is: %f " % (i, p)
            #print "x4 value is" , x4
            #print "value of p", p
            #print "value of p_new", p_new
            
            if p_new < p:
                p_new = p
                #print "new value of p_new is:  ", p_new
                x+=1
                #print "yes"
                used = True
                #print used
            else:
                #print "want to see if use is true or false  :", used
                #print "no"
              
                p_best = p_new
                if used:
                    x_best = x-1
                    #print 'best accuracy for seed %d is %f' %(i, p_best)
                    #print "best parameter so far: %d" %x_best
                break 
    #print "best accuracy for seed is %f" % p_best
    #print "best parameter so far: %d" %x_best
    return x_best, p_best


def seeding3(f1, f2, x_seed = 1, iteration = 0, x2=1, x4=1, x5=1, x6=1, step = 1, first_time = True):
   
    for i in range(1,11,step):
            
        #print "seed: ", i
        
        if iteration%2!=0:
            x= i
            x7 = x_seed
            
            #print "x4 is parameter"
            
        else:
            x = i
            x8= x_seed
            #print "x2 is parameter"
            
        p = 0.0     
        
        if first_time:
            p_new= -1
            p_best = -1
        else:
            p_new = p_best
        first_time = False
        
        used = False
        
        while True:
            if iteration%2!=0:
                x8= x
            else:
                x7 = x
            
            #print "count:  ", count
            #print "value of x4: ", x4
            p = parameters(f1, f2, 0, x2, 0, x4, x5, x6, x7, x8)
            #print "accuracy for seed %d is: %f " % (i, p)
            #print "x4 value is" , x4
            #print "value of p", p
            #print "value of p_new", p_new
            
            if p_new < p:
                p_new = p
                #print "new value of p_new is:  ", p_new
                x+=1
                #print "yes"
                used = True
                #print used
            else:
                #print "want to see if use is true or false  :", used
                #print "no"
              
                p_best = p_new
                if used:
                    x_best = x-1
                    #print 'best accuracy for seed %d is %f' %(i, p_best)
                    #print "best parameter so far: %d" %x_best
                break 
    #print "best accuracy for seed is %f" % p_best
    #print "best parameter so far: %d" %x_best
    return x_best, p_best


def sub_tuning1(f1, f2, x5=1, x6=1, x7=5, x8=1):
    x5=1; x6=1; x7=1;
    initial = True
    iteration = 0
    while True:
        
        if initial:        
            best_accuracy = -1
            p_best = -1
            x = 1
        else:
            if best_accuracy < p_best:
                best_accuracy = p_best
                x = x_best 
            else:
                print "Best accuracy  %f" % best_accuracy
                return x2_best, x4_best
                break
        iteration += 1 
        initial = False        
        x_best, p_best = seeding1(f1, f2, x, iteration)
        if iteration%2!=0:
            x4_best= x_best
            print "x4 best so far:", x4_best
        else:
            x2_best = x_best
            print "x2 best so far:", x2_best
            



def sub_tuning2(f1, f2, x2=1, x4=1, x7=5, x8=1):
    initial = True
    iteration = 0
    while True:
        
        if initial:        
            best_accuracy = -1
            p_best = -1
            x = 0
        else:
            if best_accuracy < p_best:
                best_accuracy = p_best
                x = x_best 
            else:
                print "Best accuracy  %f" % best_accuracy
                return x5_best, x6_best
                break
        iteration += 1 
        initial = False        
        x_best, p_best = seeding2(f1, f2, x, iteration, x2, x4)
        if iteration%2!=0:
            x6_best= x_best
            print "x6 best so far:", x6_best
        else:
            x5_best = x_best
            print "x5 best so far:", x5_best
    print "Best accuracy  %f" % best_accuracy

def sub_tuning3(f1, f2, x2=1, x4=1, x5=1, x6=1):
    
    initial = True
    iteration = 0
    while True:
        
        if initial:        
            best_accuracy = -1
            p_best = -1
            x = 1
        else:
            if best_accuracy < p_best:
                best_accuracy = p_best
                x = x_best 
            else:
                print "Best accuracy  %f" % best_accuracy
                return x7_best, x8_best, best_accuracy
                break
        iteration += 1 
        initial = False        
        x_best, p_best = seeding3(f1, f2, x, iteration, x2, x4, x5, x6)
        if iteration%2!=0:
            x8_best= x_best
            print "x8 best so far:", x8_best
        else:
            x7_best = x_best
            print "x7 best so far:", x7_best
            


def tune(f1, f2):
    x2, x4 = sub_tuning1(f1, f2)
    x5, x6 = sub_tuning2(f1, f2,  x2, x4)
    x7, x8, b = sub_tuning3(f1, f2, x2, x4, x5, x6)
    return 1, x2, 1, x4, x5, x6, x7, x8, b
print tune("amp.fasta", "non_amp.fasta" )


fh.write(tune("amp.fasta", "non_amp.fasta" ))
fh.close()





