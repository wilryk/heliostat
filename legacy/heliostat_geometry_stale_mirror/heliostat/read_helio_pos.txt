#%%
import pandas as pd
from matplotlib import pyplot as plt
import cv2
from scipy.spatial import distance_matrix
import numpy as np

def select_diverse_points(points, n_select):
    """
    Select diverse points using a greedy farthest-point sampling approach.
    """
    if len(points) <= n_select:
        return list(range(len(points)))
    
    selected = [np.random.randint(len(points))]
    
    for _ in range(n_select - 1):
        # Calculate minimum distance to already selected points
        dists = distance_matrix(points, points[selected])
        min_dists = dists.min(axis=1)
        
        # Select point farthest from any selected point
        next_point = np.argmax(min_dists)
        selected.append(next_point)
    
    return selected

def visualize_results(x,y,xs,ys):
    """
    Visualize detected circles on the original image
    """

    plt.figure(figsize=(12, 12))
    plt.rc('font', size=24)  
    plt.scatter(x,y,s=10)
    plt.scatter(xs,ys,s=50)
    plt.title('Heliostat positions')
    plt.axis('off')
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig('detected_circles.png', dpi=150, bbox_inches='tight')
    plt.show()

file = r"C:\Users\eadsr\Downloads\Beam down 645 field x,y centers.xlsx"

# Read the Excel file
df = pd.read_excel(file)

# Extract columns into variables
xpos = df["x (m)"].to_numpy()
ypos = df["y (m)"].to_numpy()

points = np.column_stack((xpos,ypos))
n_select = 25

selected = select_diverse_points(points, n_select)

xs = points[selected,0]
ys = points[selected,1]

visualize_results(xpos,ypos,xs,ys)

np.savetxt('positions_selected.csv', np.column_stack((xs,ys)), 
               delimiter=',', header='X (m),Y (m)', comments='')
print("\nCoordinates saved to 'positions_selected.csv'")