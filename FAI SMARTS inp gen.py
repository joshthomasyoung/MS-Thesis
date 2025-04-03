import pandas as pd
import os

# Define file path
file_path = "Data/Input Data/AK_FAIRBANKS-IAP_702610_TMY3.csv"

# Read the data with a specified encoding and only load rows 2880 to 5784
data = pd.read_csv(file_path, skiprows=18, encoding="ISO-8859-1")

# Select and transform necessary columns for rows 2880 to 5784
data = data.iloc[2880:5785, [0, 1, 6, 21, 25, 29]]  # Rows 2880-5784 (1-based) map to indices 2880-5785 (0-based)
data.columns = ["Date", "Hour", "Site Pressure (mb)", "Visibility (km)", "Precipitable Water (cm)", "Albedo"]
data["Site Pressure (mb)"] = data["Site Pressure (mb)"] / 100
data["Precipitable Water (cm)"] = data["Precipitable Water (cm)"] / 10

# Set the upper limit for precipitable water in SMARTS model
data["Precipitable Water (cm)"] = data["Precipitable Water (cm)"].clip(upper=12)

# Ensure the Date column is in datetime format and extract Year, Month, Day
data['Date'] = pd.to_datetime(data['Date'])
data['Year'], data['Month'], data['Day'] = data['Date'].dt.year, data['Date'].dt.month, data['Date'].dt.day
data["Hour"] = data["Hour"].astype(str).str.split(":").str[0].astype(int)

# Directory to save the files
output_dir = "../SMARTS_295_PC/Fairbanks INP"
os.makedirs(output_dir, exist_ok=True)  # Create the directory if it doesn't exist

# Loop through each row in the filtered data to generate files from 2880 to 5784
for i, row in enumerate(data.iterrows(), start=2880):
    row_data = row[1]  # Access the Series data for the current row

    # Create the content for each file based on the template, using the correct file number
    content = f"""'{i}'
1
{row_data['Site Pressure (mb)']:.3f} .136 .012192
1
'SAS'
0
{row_data['Precipitable Water (cm)']}
1
0
3
420
0
'S&F_RURAL'
4
{row_data['Visibility (km)']}
-1
{row_data['Albedo']}
0
300 4000 1 1366.1
2
300 4000 10
1
4
0
0
0
0
3
{row_data['Year']} {row_data['Month']} {row_data['Day']} {row_data['Hour']} 64.8155 -147.8595 -9"""

    # Define the filename and save path
    filename = f"{str(i).zfill(4)}.txt"
    file_path = os.path.join(output_dir, filename)

    # Write the content to the file
    with open(file_path, "w") as file:
        file.write(content)

print("Text files generated successfully in the Data/SMARTS INPUT directory for data points 2880 to 5784.")
