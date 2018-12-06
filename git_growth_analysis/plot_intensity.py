import scipy
import numpy as np 
import matplotlib.pyplot as plt
import os
import skimage
import time
import operator
import tkinter as tk

from statistics import mean
from scipy import optimize
from scipy.optimize import leastsq
from scipy.optimize import fmin
from skimage import io 
from numpy import array

FILE = "Plasma Variation"
MASK = 50              #150    #165
START = 25             #200   #250. 500     
END = 800
PERIOD = 50             #45    #60, 40           

# FILE = "GaN0086"
# MASK = 165              #150    #165
# START = 250             #200   #250. 500     
# PERIOD = 40             #45    #60, 40
past = time.clock()


'''
Capture the interest points that we need
input: the image and the minimum intensity value
output: the image with only dots we need
'''
def capture(semi, m_val):
        green = semi[:,:,1]             #filter out all colors except green
        mask = green > m_val
        interest = green * mask         #apply mask
        return interest

'''
Evaluate the user input mask by displaying the collected dots on the image with max and min
input: all available image intensity
output: none
'''
def evaluate_mask(intens_y):
        # intens_y = intens_y[START:]
        max_intens = max(intens_y)
        max_index = intens_y.index(max_intens)
        semi = io.imread(files[max_index])
        max_pic = capture(semi, MASK)           #gets the max pic and filter it with our current MASK
        # print(mean(max_pic))

        med_intens = np.median(intens_y)
        med_intens = intens_y.index(med_intens)
        semi = io.imread(files[med_intens])
        med_pic = capture(semi, MASK)

        min_intens = min(intens_y)
        min_index = intens_y.index(min_intens)  #obtain min index
        semi = io.imread(files[min_index])
        min_pic = capture(semi, MASK)
        
        plt.subplot(311)                      #check our MASK
        plt.imshow(max_pic)
        plt.subplot(312)
        plt.imshow(med_pic)
        plt.subplot(313)
        plt.imshow(min_pic)
        plt.show()

'''
Calculate the period and frequency of the graph by finding minimum groups 
input: mean, std, guess period, and dictionary for the data
output: calculated period
'''
def minimum_grouping(m, std, period, dictionary):
        lowest = m - 1.5*std        #this only gives us the last 5% of the intensity
        gap = period * .33
        low_x = []
        low_y = []
        for key in sorted(dictionary.keys()):
                if(dictionary[key] < lowest):
                        low_x.append(key)
                        low_y.append(dictionary[key])

        #obtain the groups for each low points
        group = []
        temp = []
        temp.append(low_x[0])
        for i in range(1, len(low_x)):
                if(low_x[i] - low_x[i-1] < gap):
                        temp.append(low_x[i])
                else:
                        # print(temp)
                        group.append(mean(temp))
                        temp.clear()
                        temp.append(low_x[i])
        group.append(mean(temp))        #this has the mean value of each group

        #find the difference of each group and return the mean of the difference
        return mean([j-i for i, j in zip(group[:-1], group[1:])])

'''
slice the data into groups of guess_period*1.3 and search for the max data in that group
input: guess_period, after start x and y data
output: the mean of the difference between maximas = mean of period
'''
def maximum_selecting(guess_period, x, y):
        i = 0
        count = 0
        temp_x = []
        temp_y = []
        x_val = []
        while (i < len(x)):
                if(count == int(guess_period*1.3)):    #guess*1.3 is when we select max
                        diction = dict(zip(temp_x, temp_y))
                        x_val.append(max(diction.items(), key=operator.itemgetter(1))[0]) #append the x val of the max y val on x_val
                        temp_x.clear()
                        temp_y.clear()
                        count = 0
                        i = int(x_val[-1] - START + guess_period*0.2)  #reset iterator to the found max's x * 1.2 to avoid mutiple counts
                else:
                        temp_x.append(x[i])
                        temp_y.append(y[i])
                        count += 1
                        i+=1
        return mean([j-i for i, j in zip(x_val[:-1], x_val[1:])])   #return the mean of the difference

'''
obtain the max 5% of the data and min 5% of the data to perform scipy fit on it
input: mean, standard deviation, and x, y data
output: numpy arrays of x coordinates and y coordinates
'''
def maximum_minimum_scipy(m, std, dictionary):
        lowest = m - 1.3*std        #this only gives us the last 5% of the intensity
        lowest_cap = m - 3*std
        low_x = []
        low_y = []
        for key in sorted(dictionary.keys()):
                if(dictionary[key] < lowest and dictionary[key] > lowest_cap):
                        low_x.append(key)
                        low_y.append(dictionary[key])
        # plt.subplot(211)
        # plt.plot(low_x, low_y, 'ro')

        highest = m + 1.3*std        #this only gives us the top 5% of the intensity
        highest_cap = m + 3*std
        high_x = []
        high_y = []
        for key in sorted(dictionary.keys()):
                if(dictionary[key] > highest and dictionary[key] < highest_cap):
                        high_x.append(key)
                        high_y.append(dictionary[key])
        # plt.plot(high_x, high_y, 'bo')
        # plt.show()
        return (array(low_x+high_x), array(low_y+high_y))

#finding the correct path
path = os.getcwd() + "/" + FILE                                                                  #this requires user input!
os.chdir(path)

#put all .txt files in files
files = []              #contains all of the .tiff files
for f in os.listdir():
    if f.endswith(".tiff"):
        files.append(f)     #put every image in files

#put all images intensity in a list, y axis
intens_y = []
check_mask = []
for f in files:
    semi = io.imread(f)
    interest = capture(semi, MASK)  #imply mask on each images                                  #this requires user input
    intens_y.append(interest.mean() * 100)  #put the average intensity in a list

#check if our mask is good
# evaluate_mask(intens_y)
/Users/joseph/Desktop/GaN/period_os
#start getting our x axis
time_x = []
for f in files:
    time_x.append(float(f[36] + f[37] + f[38] + f[39]))

plt.plot(time_x, intens_y, 'ro')             #check whole graph (sinc)
plt.show()


#----------convert into dictionary to only get the sin part of the function
diction = dict(zip(time_x, intens_y))

#sort the dic by time and only capture the parts that we are interested, aka when we start growth
x_plot = []
y_plot = []
for key in sorted(diction.keys()):
        if(key > START and key < END):                                                                        #this requires user input
                x_plot.append(key)
                y_plot.append(diction[key])
x_np = array(x_plot)
y_np = array(y_plot)
dictionary = dict(zip(x_np, y_np))

plt.plot(x_np, y_np, 'ro')            #check the START (sin)
plt.show()


#--------Guess sin fitting curve values!
the_mean = np.mean(y_np)
the_std = np.std(y_np)
guess_phase = 0
guess_period = PERIOD                                                                                 #this requires user's input
guess_freq = 1/guess_period
guess_amp = the_std*2
guess_slope = 0
# print("Guess: ", PERIOD)
data_guess = (guess_amp)*np.cos(x_np*guess_freq*2*np.pi) + the_mean     #equations for plot


#minimum grouping
#get the calculated frequency and period
min_group_period = minimum_grouping(the_mean, the_std, guess_period, dictionary)
min_group_freq = 1/min_group_period
print("Group Min:", min_group_period)
data_min_group = (guess_amp)*np.cos(x_np*min_group_freq*2*np.pi) + the_mean     #equations for plot


#maximum selecting
max_sel_period = maximum_selecting(guess_period, x_np, y_np)
max_sel_freq = 1/max_sel_period
print("Select Max: ", max_sel_period)
data_max_sel = (guess_amp)*np.cos(x_np*max_sel_freq*2*np.pi) + the_mean     #equations for plot


#perform algorithm fit!
# optimize_func = lambda x: x[0]*np.sin(x[1]*x_np+x[2]) + x[3] - y_np
optimize_func = lambda x: x[0]*np.sin(x[1]*x_np+x[2]) + the_mean + x_np*x[3]
# print(leastsq(optimize_func, [guess_amp, guess_freq, guess_phase, the_mean]))
est_amp, est_freq, est_phase, est_mean, est_slope = leastsq(optimize_func, [guess_amp, guess_freq, guess_phase, the_mean, guess_slope])[0]
fine_t = np.arange(START,max(x_np),0.1)
print("Scipy: ", 1/est_freq)
data_fit = est_amp*np.sin((fine_t*est_freq+est_phase)*2*np.pi) + the_mean + est_slope*fine_t     #equations for plot


#take max std and min std part of the graph and scipy it
mms_x, mms_y = maximum_minimum_scipy(the_mean, the_std, dictionary)
mms_func = lambda x: x[0]*np.sin(x[1]*mms_x+x[2]) + the_mean + mms_x*x[3]
mms_amp, mms_freq, mms_phase, mms_mean, est_slope = leastsq(mms_func, [guess_amp, guess_freq, guess_phase, the_mean, guess_slope])[0]
print("Max Min Scipy: ", 1/mms_freq)
data_mms = guess_amp*np.sin((fine_t*mms_freq+mms_phase)*2*np.pi) + the_mean + est_slope*fine_t    #equations for plot


# # # plot our data for user guess
plt.subplot(511)                   #for guess
plt.title('Guess')
plt.plot(x_plot, y_plot, 'ro')
plt.plot(x_plot, data_guess, color = 'm')

plt.subplot(512)                   #for calculated
plt.title('Min_Grouping')
plt.plot(x_plot, y_plot, 'ro')
plt.plot(x_plot, data_min_group, color = 'b')

plt.subplot(513)                   #for calculated
plt.title('Selecting Max')
plt.plot(x_plot, y_plot, 'ro')
plt.plot(x_plot, data_max_sel, color = 'b')

plt.subplot(514)                   #for estimate
plt.title('Scipy')
plt.plot(x_plot, y_plot, 'ro')
plt.plot(fine_t, data_fit)

plt.subplot(515)                   #for max min scy, SOLO
plt.title('Max Min Scipy')
plt.plot(x_plot, y_plot, 'ro')
plt.plot(fine_t, data_mms)
plt.show()


# -----printing out the time it takes
now = time.clock()
dif = now - past
print('time: ' + str(dif))

