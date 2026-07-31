#%%
from axicon_angle_bounds import find_valid_angles

results = find_valid_angles(h=10, x_m=5, y_m=3, x_min=2, x_max=8, plot=False)

angle_min = results['valid_ranges'][0][0]
angle_max = results['valid_ranges'][0][1]
print(f"Minimum angle is {angle_min} rad, and maximum angle is {angle_max}")
