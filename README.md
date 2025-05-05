This project contains the important files for my MS Thesis and Python model. The Python model is structured to be adaptable for different locations. The files are described below:

Young_Thesis_April24-V7.pdf - This is my MS Thesis. This thesis fully explains the development of the Python model including methodology, model improvements, results, and discussion. This thesis also explores an application of the model used to optimize HVAC systems in Fairbanks, AK.

Fairbanks Model Home V4.py - This is the main file for the Python model.

FAI SMARTS inp gen.py - This Python file is used to generate a year's worth of input files for SMARTS_295. Note that the use of SMARTS_295 is required with this model.

SSI Generator.bat - This windows batch script is used to iteratively run SMARTS_295 application on desktop, 8760 times corresponding to every hour for one year, using each of the input files generated from FAI SMARTS inp gen.py and generating/storing an output file for each. The SSI data must be generated and stored correctly before using the main file. The main python file pulls directly from this generated data.
