#!/usr/bin/python

"""
    Script to automatically measure the fluorescence intensity of vagus axon branches.  Input is CSV files from Plot Profile function in ImageJ.
    ***NOTE: IF YOU ENCOUNTER AN "ASCII" OR "UNICODE" ERROR, CSV FILES FROM IMAGEJ MUST BE OPENED AND SAVED IN CSV UTF-8 FORMAT, WHICH WILL SOLVE THE PROBLEM
    ***NOTE: CSV FILES MUST BE MOVED OFF OF THE SERVER AND ONTO THE COMPUTER HARD DRIVE, BEFORE RUNNING THIS SCRIPT.
    ***NOTE: CSV FILES MUST BE NAMED "...green.csv" and "...red.csv"
    
"""

import sys
import csv
import os
import statistics
import pandas as pd
import matplotlib.pyplot as plt


#define source folder *****THIS NEEDS TO BE CHANGED FOR EACH EXPERIMENT*****
#source = "/Volumes/aisabell/moenslab/Adam/Confocal/2018-2-23 isl1Kaede crizotinib 3dpf/Analysis 2018-02-23/"
source = "/Users/aisabell/Documents/Moens Local/Analysis/2019-07-23 dpysl3 taeer/"

#generate output csv file, write header row
summary_file = source + "python_results.csv"
with open(summary_file, 'w') as t:
    writer = csv.writer(t)
    writer.writerow(('file', 'number of peaks', 'green background', 'red background', 'signal', 'branch number', 'position range', 'number of positions', 'green sum intensity', 'green max intensity', 'green average intensity', 'red sum intensity', 'red max intensity', 'red average intensity'))

#command to run OP analysis on all .csv files within the source folder
for filename in os.listdir(source):
    if filename.endswith("green.csv") and "python" not in filename:
        green_filename = source + filename
        red_filename = source + filename.replace("green", "red")
        short_filename = green_filename.replace("green.csv", "")
        print(short_filename)
        
        #read in source file, write each column to a list (there will be two columns: column 1 is the position, column 2 is the fluorescence value.)        
        green_df = pd.read_csv(green_filename)
        pos_list = green_df['X'].tolist()
        green_fluor_list = green_df['Y'].tolist()
        red_df = pd.read_csv(red_filename)
        red_fluor_list = red_df['Y'].tolist()      

        #values are read into lists as strings.  convert to floats.
        green_fluor_list = list(map(float, green_fluor_list))
        red_fluor_list = list(map(float, red_fluor_list))

        #measure the mean of the first 10 values in fluor_list to get the background level.  set level of real signal to 3 times background.
        background = statistics.mean(green_fluor_list[0:10])
        red_background = statistics.mean(red_fluor_list[0:10])
        signal = background * 3
        if signal < 1:
            signal = 1

        #iterate over fluor_list, identifying groups of values that represent real signal peaks.  pass each of these groups to a list, and when the grouping is complete, save that data, for both fluor_list and pos_list.
        green_val_fluor_list = []
        val_pos_list = []
        green_val_fluor_list_master = []
        val_pos_list_master = []
        red_val_fluor_list = []
        red_val_fluor_list_master = []
        val = 0
        while val < len(green_fluor_list):
            if green_fluor_list[val] > signal:
                green_val_fluor_list.append(green_fluor_list[val])
                val_pos_list.append(pos_list[val])
                red_val_fluor_list.append(red_fluor_list[val])                
            val = val + 1
            try:
                if green_fluor_list[val] < signal:
                    if len(green_val_fluor_list) >= 15:
                        val_pos_list_master.append(val_pos_list)
                        green_val_fluor_list_master.append(green_val_fluor_list)
                        red_val_fluor_list_master.append(red_val_fluor_list)
                    green_val_fluor_list = []
                    red_val_fluor_list = []
                    val_pos_list = []
            except IndexError:
                if len(green_val_fluor_list) >=15:
                    val_pos_list_master.append(val_pos_list)
                    green_val_fluor_list_master.append(green_val_fluor_list)
                    red_val_fluor_list_master.append(red_val_fluor_list)    

        #now we've got master lists, which contain lists of fluor and pos values corresponding to peaks of intensity.  here we take each of those peaks, analyze it, and write the data to our output csv.
        peak = 0
        maxvals = []
        while peak < len(green_val_fluor_list_master):
            maxvals.append(max(green_val_fluor_list_master[peak]))
            num_peaks = len(green_val_fluor_list_master)
            branch_number = peak + 4
            peak_length = len(green_val_fluor_list_master[peak]) - 1
            num_pos = peak_length + 1
            pos_range = str(val_pos_list_master[peak][0]) + " : " + str(val_pos_list_master[peak][peak_length])

            green_sum_intensity = sum(green_val_fluor_list_master[peak])
            green_max_intensity = max(green_val_fluor_list_master[peak])
            green_mean_intensity =sum(green_val_fluor_list_master[peak]) / (peak_length + 1) 
            red_sum_intensity = sum(red_val_fluor_list_master[peak])
            red_max_intensity = max(red_val_fluor_list_master[peak])
            red_mean_intensity =sum(red_val_fluor_list_master[peak]) / (peak_length + 1) 
            with open(summary_file, 'a') as t:
                writer = csv.writer(t)
                writer.writerow((short_filename, num_peaks, background, red_background, signal, branch_number, pos_range, num_pos, green_sum_intensity, green_max_intensity, green_mean_intensity, red_sum_intensity, red_max_intensity, red_mean_intensity))
            peak = peak + 1

        #now we want to generate an output graph showing the green peaks and the regions considered peaks by our software.
        #add column Z to df, populate with zeros
        dflen = len(green_df['X'])
        z = [0] * dflen
        green_df['Z'] = z
        #for each position counted in a peak, change value of Z in df to max value fluorescent intensities.
        maxvalue = max(maxvals)
        for i in green_val_fluor_list_master:
            for t in i:
                green_df.loc[green_df.Y == t, 'Z'] = maxvalue/3
        #graph real peaks and called peaks.
        green_df.columns = ['X', 'Green Fluorescence', 'Called Peaks']
        green_df.plot(x = 'X', y = ['Green Fluorescence', 'Called Peaks'])
        plt.savefig(short_filename + " called_peaks.jpg")
        plt.close()
        print("saved")