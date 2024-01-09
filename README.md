MATLAB Data Analysis Code for EMG Research

This repository contains MATLAB scripts developed for the analysis of Electromyography (EMG) data as part of my third-year project at the University of Oxford. The scripts are designed to process and analyze EMG data, focusing on segment extraction and optimal window length determination.

Files in the Repository

experimental_data.m
Description: This is the main script for the project. It processes EMG data and performs initial analyses.
Usage: Run this script in MATLAB to start the data analysis process. Ensure that your EMG data files are in the appropriate format and location as expected by the script.
extract_segment.m
Description: Contains functions that extract segments of EMG data. These functions allow for the specification of window length, start time, and sensing frequency.
Usage: This script is called by experimental_data.m. It can also be used independently for specific data segment extraction tasks.
window_length.m
Description: A separate script for analyzing the optimal window length for EMG data analysis.
Usage: Run this script to determine the most effective window length for your EMG data analysis. This can help in optimizing the data processing in experimental_data.m.
Getting Started

To use these scripts, clone this repository to your local machine. Ensure you have MATLAB installed and that your EMG data files are properly formatted and located in the directory expected by the scripts.

Prerequisites

MATLAB (version [specify version if necessary])
EMG data files in [specify expected format, e.g., .csv, .mat]
Running the Scripts

Open MATLAB.
Navigate to the directory containing the cloned repository.
Run experimental_data.m for a full analysis.
Optionally, run window_length.m to adjust the window length based on your specific dataset.

Authors

Jiaqi Li - Initial work 
