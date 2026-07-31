import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

def calculate_gamma(a, h, x_m, y_m):
    """
    Calculate gamma given angle a and parameters.
    
    gamma = (m*x_m - y_m) / (m - tan(a))
    where m = (h*cos(2*a) - y_m) / (-h*sin(2*a) - x_m)
    
    Parameters:
    -----------
    a : float
        Angle in radians
    h, x_m, y_m : float
        Problem parameters
        
    Returns:
    --------
    float : gamma value
    """
    m = (h*np.cos(2*a) - y_m) / (-h*np.sin(2*a) - x_m)
    gamma = (m*x_m - y_m) / (m - np.tan(a))
    return gamma


def find_valid_angles(h, x_m, y_m, x_min, x_max, 
                      a_range=(0.01, np.pi/2 - 0.01),
                      n_points=2000,
                      plot=False):
    """
    Find valid angle ranges where x_min < gamma < x_max.
    
    Parameters:
    -----------
    h, x_m, y_m : float
        Problem parameters
    x_min, x_max : float
        Bounds for gamma
    a_range : tuple, optional
        (min_angle, max_angle) in radians to search
    n_points : int, optional
        Number of points for numerical evaluation
    plot : bool, optional
        If True, generate visualization plots
        
    Returns:
    --------
    dict : Dictionary containing:
        - 'valid_ranges': list of tuples [(a_start, a_end), ...] in radians
        - 'valid_ranges_deg': same ranges in degrees
        - 'boundaries': list of boundary angles in radians
        - 'boundaries_deg': list of boundary angles in degrees
    """
    
    # Sample angles
    angles = np.linspace(a_range[0], a_range[1], n_points)
    gamma_values = []
    valid_angles = []
    
    # Calculate gamma for all angles
    for a in angles:
        try:
            gamma = calculate_gamma(a, h, x_m, y_m)
            if np.isfinite(gamma):
                gamma_values.append(gamma)
                valid_angles.append(a)
            else:
                gamma_values.append(np.nan)
                valid_angles.append(a)
        except:
            gamma_values.append(np.nan)
            valid_angles.append(a)
    
    gamma_array = np.array(gamma_values)
    angles_array = np.array(valid_angles)
    
    # Find where constraint is satisfied: x_min < gamma < x_max
    satisfied = (gamma_array > x_min) & (gamma_array < x_max)
    
    # Find boundaries (where gamma crosses x_min or x_max)
    boundaries = []
    
    # Check for crossings with x_min
    diff_min = gamma_array - x_min
    for i in range(len(diff_min)-1):
        if np.isfinite(diff_min[i]) and np.isfinite(diff_min[i+1]):
            if diff_min[i] * diff_min[i+1] < 0:  # Sign change
                try:
                    a_boundary = brentq(lambda a: calculate_gamma(a, h, x_m, y_m) - x_min,
                                       angles_array[i], angles_array[i+1])
                    boundaries.append(('x_min', a_boundary))
                except:
                    pass
    
    # Check for crossings with x_max
    diff_max = gamma_array - x_max
    for i in range(len(diff_max)-1):
        if np.isfinite(diff_max[i]) and np.isfinite(diff_max[i+1]):
            if diff_max[i] * diff_max[i+1] < 0:  # Sign change
                try:
                    a_boundary = brentq(lambda a: calculate_gamma(a, h, x_m, y_m) - x_max,
                                       angles_array[i], angles_array[i+1])
                    boundaries.append(('x_max', a_boundary))
                except:
                    pass
    
    # Sort boundaries by angle
    boundaries.sort(key=lambda x: x[1])
    
    # Determine valid ranges
    valid_ranges = []
    
    if len(boundaries) == 0:
        # No boundaries - check if entire range is valid
        mid_idx = len(angles_array) // 2
        if satisfied[mid_idx]:
            valid_ranges = [(a_range[0], a_range[1])]
    else:
        # Test regions between boundaries
        test_angles = [a_range[0]]
        for i in range(len(boundaries)-1):
            test_angles.append((boundaries[i][1] + boundaries[i+1][1]) / 2)
        test_angles.append(a_range[1])
        
        for i in range(len(test_angles)):
            try:
                gamma_test = calculate_gamma(test_angles[i], h, x_m, y_m)
                if x_min < gamma_test < x_max:
                    # This region is valid, find its bounds
                    if i == 0:
                        start = a_range[0]
                    else:
                        start = boundaries[i-1][1]
                    
                    if i == len(test_angles) - 1:
                        end = a_range[1]
                    else:
                        end = boundaries[i][1]
                    
                    valid_ranges.append((start, end))
            except:
                pass
    
    # Plotting
    if plot:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
        
        # Plot 1: Gamma vs angle
        ax1.plot(np.degrees(angles_array), gamma_values, 'b-', linewidth=2, label='γ(a)')
        ax1.axhline(y=x_min, color='r', linestyle='--', linewidth=2, label=f'x_min = {x_min}')
        ax1.axhline(y=x_max, color='r', linestyle='--', linewidth=2, label=f'x_max = {x_max}')
        
        # Highlight valid regions
        ax1.fill_between(np.degrees(angles_array), x_min, x_max,
                         where=satisfied, alpha=0.3, color='green',
                         label='Valid region')
        
        # Mark boundaries
        for boundary_type, a_bound in boundaries:
            gamma_bound = calculate_gamma(a_bound, h, x_m, y_m)
            ax1.plot(np.degrees(a_bound), gamma_bound, 'ro', markersize=8)
        
        ax1.set_xlabel('Angle a (degrees)', fontsize=12)
        ax1.set_ylabel('γ value', fontsize=12)
        ax1.set_title(f'γ vs Angle (h={h}, x_m={x_m}, y_m={y_m})', 
                     fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=10)
        
        # Set reasonable y-limits
        valid_gamma = gamma_array[np.isfinite(gamma_array)]
        if len(valid_gamma) > 0:
            y_center = (x_min + x_max) / 2
            y_range = max(x_max - x_min, np.percentile(valid_gamma, 75) - np.percentile(valid_gamma, 25))
            ax1.set_ylim([y_center - 2*y_range, y_center + 2*y_range])
        
        # Plot 2: Constraint satisfaction
        constraint_diff_min = gamma_array - x_min
        constraint_diff_max = x_max - gamma_array
        
        ax2.plot(np.degrees(angles_array), constraint_diff_min, 'b-', 
                linewidth=1.5, label='γ - x_min', alpha=0.7)
        ax2.plot(np.degrees(angles_array), constraint_diff_max, 'r-', 
                linewidth=1.5, label='x_max - γ', alpha=0.7)
        ax2.axhline(y=0, color='k', linestyle='-', linewidth=1)
        
        ax2.fill_between(np.degrees(angles_array), 0, 
                         np.minimum(constraint_diff_min, constraint_diff_max),
                         where=satisfied, alpha=0.3, color='green',
                         label='Both constraints satisfied')
        
        ax2.set_xlabel('Angle a (degrees)', fontsize=12)
        ax2.set_ylabel('Constraint margin', fontsize=12)
        ax2.set_title('Constraint Satisfaction (both must be > 0)', 
                     fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=10)
        
        plt.tight_layout()
        plt.show()
    
    # Prepare results
    results = {
        'valid_ranges': valid_ranges,
        'valid_ranges_deg': [(np.degrees(start), np.degrees(end)) 
                             for start, end in valid_ranges],
        'boundaries': [b[1] for b in boundaries],
        'boundaries_deg': [np.degrees(b[1]) for b in boundaries],
        'boundary_types': [b[0] for b in boundaries]
    }
    
    return results


# =============================================================================
# EXAMPLE USAGE
# =============================================================================
if __name__ == "__main__":
    # Define parameters
    h = 6.0
    x_m = 40
    y_m = -20
    x_min = 1.0
    x_max = 8.0
    
    print("="*70)
    print("FINDING VALID ANGLE RANGES")
    print("="*70)
    print(f"Parameters: h={h}, x_m={x_m}, y_m={y_m}")
    print(f"Constraint: {x_min} < γ < {x_max}")
    print(f"where γ = (m*x_m - y_m)/(m - tan(a))")
    print(f"and m = (h*cos(2*a) - y_m)/(-h*sin(2*a) - x_m)")
    print("="*70)
    
    # Find valid angles with plotting
    results = find_valid_angles(h, x_m, y_m, x_min, x_max, plot=True)
    
    # Display results
    print("\n" + "-"*70)
    print("RESULTS")
    print("-"*70)
    
    if len(results['boundaries']) > 0:
        print(f"\nFound {len(results['boundaries'])} boundary point(s):")
        for i, (btype, a_rad, a_deg) in enumerate(zip(results['boundary_types'],
                                                       results['boundaries'],
                                                       results['boundaries_deg'])):
            print(f"  {i+1}. Crosses {btype} at a = {a_rad:.6f} rad = {a_deg:.4f}°")
    
    if len(results['valid_ranges']) > 0:
        print(f"\n✓ Found {len(results['valid_ranges'])} valid angle range(s):")
        for i, ((start, end), (start_deg, end_deg)) in enumerate(zip(results['valid_ranges'],
                                                                      results['valid_ranges_deg'])):
            print(f"\n  Range {i+1}:")
            print(f"    {start:.6f} ≤ a ≤ {end:.6f} radians")
            print(f"    {start_deg:.4f}° ≤ a ≤ {end_deg:.4f}°")
            
            # Verify by checking a point in the middle
            a_test = (start + end) / 2
            gamma_test = calculate_gamma(a_test, h, x_m, y_m)
            print(f"    Verification at midpoint: γ({np.degrees(a_test):.2f}°) = {gamma_test:.4f}")
    else:
        print("\n✗ No valid angle ranges found for the given constraints.")
    
    print("\n" + "="*70)