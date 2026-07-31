import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, brentq

def calculate_m(a, h, x_m, y_m):
    """Calculate m given angle a and parameters"""
    return (h*np.cos(2*a) - y_m) / (-h*np.sin(2*a) - x_m)

def calculate_lhs(a, h, x_m, y_m):
    """Calculate left-hand side of inequality"""
    m = calculate_m(a, h, x_m, y_m)
    return (m*x_m - y_m) / (m - np.tan(a))

def inequality_function(a, h, x_m, y_m, x_min):
    """Returns LHS - x_min (positive when inequality satisfied)"""
    return calculate_lhs(a, h, x_m, y_m) - x_min

# =============================================================================
# EXAMPLE: Solve for specific parameter values
# =============================================================================
print("="*70)
print("NUMERICAL SOLUTION FOR INEQUALITY: (m*x_m - y_m)/(m-tan(a)) > x_min")
print("where m = (h*cos(2*a) - y_m)/(-h*sin(2*a) - x_m)")
print("="*70)

# Set your parameters here
h = 10.0
x_m = 5.0
y_m = 3.0
x_min = 2.0

print(f"\nParameters: h={h}, x_m={x_m}, y_m={y_m}, x_min={x_min}")

# =============================================================================
# STEP 1: Visualize the inequality across angles
# =============================================================================
print("\n" + "-"*70)
print("STEP 1: Visualizing the inequality")
print("-"*70)

angles = np.linspace(0.01, np.pi/2 - 0.01, 1000)
lhs_values = []
valid_angles = []

for a in angles:
    try:
        lhs = calculate_lhs(a, h, x_m, y_m)
        if np.isfinite(lhs):  # Check for valid values
            lhs_values.append(lhs)
            valid_angles.append(a)
        else:
            lhs_values.append(np.nan)
            valid_angles.append(a)
    except:
        lhs_values.append(np.nan)
        valid_angles.append(a)

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# Plot 1: LHS vs angle
ax1.plot(np.degrees(valid_angles), lhs_values, 'b-', linewidth=2, label='LHS')
ax1.axhline(y=x_min, color='r', linestyle='--', linewidth=2, label=f'x_min = {x_min}')
ax1.fill_between(np.degrees(valid_angles), x_min, lhs_values, 
                  where=np.array(lhs_values) > x_min, alpha=0.3, color='green',
                  label='Inequality satisfied')
ax1.set_xlabel('Angle a (degrees)', fontsize=12)
ax1.set_ylabel('LHS value', fontsize=12)
ax1.set_title('Left-hand side vs Angle', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=10)
ax1.set_ylim([x_min - 5, x_min + 10])

# Plot 2: Inequality difference (LHS - x_min)
diff = np.array(lhs_values) - x_min
ax2.plot(np.degrees(valid_angles), diff, 'g-', linewidth=2)
ax2.axhline(y=0, color='k', linestyle='-', linewidth=1)
ax2.fill_between(np.degrees(valid_angles), 0, diff, 
                  where=diff > 0, alpha=0.3, color='green', label='> 0 (satisfied)')
ax2.fill_between(np.degrees(valid_angles), 0, diff, 
                  where=diff < 0, alpha=0.3, color='red', label='< 0 (not satisfied)')
ax2.set_xlabel('Angle a (degrees)', fontsize=12)
ax2.set_ylabel('LHS - x_min', fontsize=12)
ax2.set_title('Inequality Difference (positive = satisfied)', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=10)

plt.tight_layout()
plt.savefig('inequality_solution.png', dpi=150, bbox_inches='tight')
print("\nPlot saved as 'inequality_solution.png'")

# =============================================================================
# STEP 2: Find boundary angles where LHS = x_min
# =============================================================================
print("\n" + "-"*70)
print("STEP 2: Finding boundary angles (where LHS = x_min)")
print("-"*70)

# Find sign changes to locate roots
diff_array = np.array(lhs_values) - x_min
sign_changes = []

for i in range(len(diff_array)-1):
    if np.isfinite(diff_array[i]) and np.isfinite(diff_array[i+1]):
        if diff_array[i] * diff_array[i+1] < 0:  # Sign change
            sign_changes.append((valid_angles[i], valid_angles[i+1]))

print(f"\nFound {len(sign_changes)} potential boundary region(s)")

boundaries = []
for i, (a_low, a_high) in enumerate(sign_changes):
    try:
        # Use Brent's method for robust root finding
        a_boundary = brentq(inequality_function, a_low, a_high, 
                           args=(h, x_m, y_m, x_min))
        boundaries.append(a_boundary)
        print(f"\nBoundary {i+1}:")
        print(f"  a = {a_boundary:.6f} radians = {np.degrees(a_boundary):.4f} degrees")
        print(f"  Verification: LHS = {calculate_lhs(a_boundary, h, x_m, y_m):.6f}")
    except:
        print(f"\nBoundary {i+1}: Failed to converge")

# =============================================================================
# STEP 3: Determine solution regions
# =============================================================================
print("\n" + "-"*70)
print("STEP 3: Solution regions where inequality is satisfied")
print("-"*70)

if len(boundaries) > 0:
    # Test regions between boundaries
    test_points = [0.1]  # Start
    test_points.extend([(boundaries[i] + boundaries[i+1])/2 
                       for i in range(len(boundaries)-1)])
    if boundaries[-1] < np.pi/2 - 0.1:
        test_points.append(boundaries[-1] + 0.1)
    
    print("\nTesting regions:")
    for i, test_a in enumerate(test_points):
        if test_a < np.pi/2:
            try:
                lhs_test = calculate_lhs(test_a, h, x_m, y_m)
                satisfied = lhs_test > x_min
                print(f"\nRegion {i+1}: a ≈ {np.degrees(test_a):.2f}°")
                print(f"  LHS = {lhs_test:.4f}, x_min = {x_min}")
                print(f"  Inequality satisfied: {satisfied}")
            except:
                print(f"\nRegion {i+1}: Undefined or singularity")
else:
    print("\nNo boundaries found - inequality may be always satisfied or never satisfied")
    print("Check the plot for visual confirmation")

print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"For the given parameters, the inequality (m*x_m - y_m)/(m-tan(a)) > {x_min}")
print("is satisfied in the green regions shown in the plot.")
print("\nTo use with different parameters, modify the values at the top of the script.")

plt.show()