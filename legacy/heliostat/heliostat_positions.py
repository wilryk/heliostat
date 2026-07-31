#%%
import pandas as pd

def read_helio_positions(file):
    # Read the csv file
    df = pd.read_csv(file,header=0)

    # Extract columns into variables
    xpos = df["X (m)"].to_numpy()
    ypos = df["Y (m)"].to_numpy()


    print("\nRead in heliostat coordinates.")
    return xpos, ypos