# -*- coding: utf-8 -*-
"""
QuadoaIrradiance Class - A comprehensive class for optical ray tracing and irradiance analysis using Quadoa
Modified to preserve NaN values and support Quadoa-saved ray trace files
"""

from quadoa import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams as rcp
import os
import time
import json
import pandas as pd
from matplotlib import cm
from matplotlib.colors import ListedColormap
import warnings

rcp.update({'font.size':12})

class QuadoaIrradiance:
    def __init__(self, model_folder, model_filename, output_directory, 
                 quadoa_base_folder="C:/Program Files/Quadoa", watts_of_source=34212, surf=None):
        """
        Initialize QuadoaIrradiance class
        
        Parameters:
        -----------
        model_folder : str
            Path to folder containing the optical model
        model_filename : str
            Name of the .optx model file
        output_directory : str
            Directory for saving results
        quadoa_base_folder : str
            Path to Quadoa installation folder (default: "C:/Program Files/Quadoa")
        watts_of_source : float
            Power of the light source in watts
        surf : int or None
            Surface number for ray tracing (None = image plane)
        """
        self.quadoa_base_folder = quadoa_base_folder
        self.model_folder = model_folder
        self.model_filename = model_filename
        self.output_directory = output_directory
        self.watts_of_source = watts_of_source
        self.surf = surf
        
        # Initialize Quadoa core
        self.core = QuadoaCore()
        
        # Create output directories
        self._create_directories()
        
        # Load materials and model
        self._load_materials_and_model()
        
        # Store ray data and tracing parameters
        self.ray_data = {}
        self.ray_dir_data ={}
        self.metadata = {}
        self.nr_rays = None
        self.n_traces = None
        
    def _create_directories(self):
        """Create necessary output directories"""
        directories = [
            self.output_directory,
            os.path.join(self.output_directory, 'ray_trace'),
            os.path.join(self.output_directory, 'irrad')
        ]
        
        for directory in directories:
            if not os.path.exists(directory):
                os.makedirs(directory)
                print(f"Created directory: {directory}")
    
    def _load_materials_and_model(self):
        """Load material catalogs and optical model"""
        try:
            # Load material catalogs
            schott_path = os.path.join(self.quadoa_base_folder, "glass", "SCHOTT.glas")
            misc_path = os.path.join(self.quadoa_base_folder, "glass", "Experimental Misc.glas")
            
            if os.path.exists(schott_path):
                self.core.loadMaterialFile(schott_path)
            else:
                print(f"Warning: SCHOTT.glas not found at {schott_path}")
                
            if os.path.exists(misc_path):
                self.core.loadMaterialFile(misc_path)
            else:
                print(f"Warning: Experimental Misc.glas not found at {misc_path}")
            
            # Load model file
            model_path = os.path.join(self.model_folder, self.model_filename)
            if os.path.exists(model_path):
                self.core.loadModelFile(model_path)
                self.core.applyChangesAndInitModel()
                print(f"Successfully loaded model: {model_path}")
            else:
                print(f"Model file not found: {model_path}")
                print(f"Model folder: {self.model_folder}")
                print(f"Model filename: {self.model_filename}")
                print(f"Files in model folder:")
                try:
                    files = os.listdir(self.model_folder)
                    for file in files:
                        print(f"  - {file}")
                except Exception as e:
                    print(f"  Could not list files: {e}")
                raise FileNotFoundError(f"Model file not found: {model_path}")
                
        except Exception as e:
            print(f"Error loading materials or model: {e}")
            raise
    
    def _rot_mat(self, theta):
        """Create rotation matrix for given angle"""
        return np.array([[np.cos(theta), -np.sin(theta)],
                        [np.sin(theta),  np.cos(theta)]])
    
    def _save_metadata(self, metadata):
        """Save metadata to JSON file"""
        metadata_path = os.path.join(self.output_directory, 'metadata.json')
        
        # Load existing metadata if it exists
        existing_metadata = {}
        if os.path.exists(metadata_path):
            try:
                with open(metadata_path, 'r') as f:
                    existing_metadata = json.load(f)
            except:
                pass
        
        # Update with new metadata
        timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
        existing_metadata[timestamp] = metadata
        
        # Save updated metadata
        with open(metadata_path, 'w') as f:
            json.dump(existing_metadata, f, indent=2)
    
    def _remove_nan_values(self, ray_positions):
        """Remove NaN values from ray positions for analysis"""
        if ray_positions.size == 0:
            return ray_positions
        # Remove columns where any value is NaN
        return ray_positions[:, ~np.isnan(ray_positions).any(axis=0)]
    
    def _count_valid_rays(self, ray_positions):
        """Count number of valid (non-NaN) rays"""
        if ray_positions.size == 0:
            return 0
        return np.sum(~np.isnan(ray_positions).any(axis=0))
    
    def set_parameters(self,solaz,solze,posx,posy,rot_az,rot_el,c3,c4,c5,config_num=0):
        self.core.setMulticonfParam("solaz",config_num,solaz)
        self.core.setMulticonfParam("solze",config_num,solze)
        self.core.setMulticonfParam("posx",config_num,posx)
        self.core.setMulticonfParam("posy",config_num,posy)
        self.core.setMulticonfParam("rot_az",config_num,rot_az)
        self.core.setMulticonfParam("rot_el",config_num,rot_el)
        self.core.setMulticonfParam("c3",config_num,c3)
        self.core.setMulticonfParam("c4",config_num,c4)
        self.core.setMulticonfParam("c5",config_num,c5)


    def quick_irradiance(self, seq_nr, grid_size, nr_rays, config_num=0, plot=True, save=True):
        """
        Use built-in Quadoa irradiance function
        
        Parameters:
        -----------
        seq_nr : int
            Sequence number to trace
        grid_size : int
            Grid resolution for irradiance image
        nr_rays : int
            Number of rays to trace
        config_num : int
            Configuration number to use
        plot : bool
            Whether to generate and save a plot
        save : bool
            Whether to save irradiance data to CSV
        
        Returns:
        --------
        irrad_img : ImageObjectD
            Quadoa irradiance image object
        """
        try:
            # Set ray distribution count
            self.core.setRayDistributionCount1(seq_nr, nr_rays)
            
            self.core.setConfig(config_num)
            self.core.traceRays(seq_nr, 0, 0)
            irrad_img = self.core.getIrradianceImageIncoherent(seq_nr, 0, 0, grid_size)
            print(f"Quick irradiance completed for config {config_num} with {nr_rays:,} rays")
            
            if plot or save:
                # Extract data from Quadoa irradiance object
                xmax = irrad_img.getXMax()
                xmin = irrad_img.getXMin()
                ymax = irrad_img.getYMax()
                ymin = irrad_img.getYMin()
                irrad_data = irrad_img.getData()
                
                print(f"Irradiance bounds: X({xmin:.2f} to {xmax:.2f}), Y({ymin:.2f} to {ymax:.2f})")
                
                # Convert to numpy array and reshape if needed
                if hasattr(irrad_data, '__len__'):
                    irrad_array = np.array(irrad_data).reshape(grid_size, grid_size)
                else:
                    print("Warning: Could not extract irradiance data array")
                    return irrad_img
                
                if save:
                    self._save_quick_irradiance_data(config_num, irrad_array, 
                                                   xmin, xmax, ymin, ymax, grid_size)
                
                if plot:
                    self._plot_quick_irradiance(config_num, irrad_array, 
                                              xmin, xmax, ymin, ymax, grid_size)
                
                # Store quick irradiance map for encircled energy analysis
                if not hasattr(self, 'irradiance_maps'):
                    self.irradiance_maps = {}
                    
                self.irradiance_maps[f"quick_{config_num}"] = {
                    'irrad_map': irrad_array,
                    'irrad_scal': 1.0,  # Quick irradiance is already scaled
                    'D': max(abs(xmin), abs(xmax), abs(ymin), abs(ymax)),
                    'grid_size': grid_size,
                    'xmin': xmin,
                    'xmax': xmax,
                    'ymin': ymin,
                    'ymax': ymax
                }
            
            return irrad_img
            
        except Exception as e:
            print(f"Error in quick_irradiance: {e}")
            raise
    
    def _save_quick_irradiance_data(self, config_num, irrad_array, xmin, xmax, ymin, ymax, grid_size):
        """Save quick irradiance data to CSV"""
        # Create coordinate grids
        x_coords = np.linspace(xmin, xmax, grid_size)
        y_coords = np.linspace(ymin, ymax, grid_size)
        X, Y = np.meshgrid(x_coords, y_coords)
        
        # Flatten arrays for CSV output
        irrad_flat = irrad_array.flatten()
        x_flat = X.flatten() / 1000  # Convert mm to m
        y_flat = Y.flatten() / 1000  # Convert mm to m
        z_flat = np.zeros_like(x_flat)
        
        # Combine data
        data = np.column_stack((irrad_flat, x_flat, y_flat, z_flat))
        
        # Save to CSV
        filename = os.path.join(self.output_directory, 'irrad', f'qirrad_{config_num + 1}.csv')
        header = "Flux (W/m2),X (m),Y (m),Z (m)"
        np.savetxt(filename, data, delimiter=',', fmt='%3.8f', header=header, comments='')
        print(f"Saved quick irradiance data: {filename}")
    
    def _plot_quick_irradiance(self, config_num, irrad_array, xmin, xmax, ymin, ymax, grid_size, dpi=300):
        """Plot and save quick irradiance data"""
        # Setup colormap (same as regular irradiance plots)
        magma = cm.get_cmap('magma', 128)
        list1 = magma.colors[30:31]
        list2 = magma.colors[45:127]
        colormap = ListedColormap(np.concatenate((list1, list2), axis=0))
        
        # Create figure
        plt.rc('font',size=14)
        fig, ax = plt.subplots(figsize=(8, 8), dpi=dpi)
        
        # Convert irradiance to kW/m^2 (assuming Quadoa gives W/m^2)
        scaled_irrad = irrad_array / 1000
        
        # Plot
        img = ax.imshow(scaled_irrad, origin='lower', cmap=colormap, 
                       vmin=0, vmax=scaled_irrad.max(), 
                       extent=[xmin, xmax, ymin, ymax])
        
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)
        plt.xlabel('X (mm)')
        plt.ylabel('Y (mm)')
        plt.title(f'Quick Irradiance - Config {config_num + 1}')
        
        # Add colorbar
        cbar = fig.colorbar(img, ax=ax)
        cbar.set_label('Irradiance (kW/m²)', fontsize=12, labelpad=10)
        
        ax.set_aspect('equal', adjustable='box')
        
        # Calculate total power
        pixel_area_mm2 = ((xmax - xmin) / grid_size) * ((ymax - ymin) / grid_size)
        pixel_area_m2 = pixel_area_mm2 / 1000000
        total_power = irrad_array.sum() * pixel_area_m2 / 1000  # Convert to kW
        
        # Save plot
        filename = os.path.join(self.output_directory, 'irrad', 
                              f'qirrad_{config_num + 1}_P_{total_power:.1f}kW.png')
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
        print(f"Saved quick irradiance plot: {filename}")
        
        plt.show()
    
    def calculate_encircled_energy(self, config_keys=None, sz=101, save=True, plot=True, D=None):
        """
        Calculate and plot encircled energy from stored irradiance maps
        
        Parameters:
        -----------
        config_keys : list or None
            List of configuration keys to analyze (None = all stored maps)
        sz : int
            Number of radius points for encircled energy calculation
        save : bool
            Whether to save encircled energy data to CSV
        plot : bool
            Whether to generate and save plots
        
        Returns:
        --------
        encircled_energy_data : dict
            Dictionary containing radius and encircled energy arrays for each config
        """
        if not self.irradiance_maps:
            print("No irradiance maps found. Run generate_irradiance_plots() or quick_irradiance() first.")
            return {}
        
        # Determine which configs to process
        if config_keys is None:
            config_keys = list(self.irradiance_maps.keys())
        
        encircled_energy_data = {}
        
        for config_key in config_keys:
            if config_key not in self.irradiance_maps:
                print(f"Warning: Configuration '{config_key}' not found in stored irradiance maps")
                continue
            
            print(f"Calculating encircled energy for {config_key}")
            
            # Get stored irradiance data
            irrad_data = self.irradiance_maps[config_key]
            irrad_map = irrad_data['irrad_map']
            irrad_scal = irrad_data['irrad_scal']
            grid_size = irrad_data['grid_size']
            
            # Handle different map types (quick vs custom)
            if 'xmin' in irrad_data:  # Quick irradiance
                xmin, xmax = irrad_data['xmin'], irrad_data['xmax']
                ymin, ymax = irrad_data['ymin'], irrad_data['ymax']
                w = max(xmax - xmin, ymax - ymin)  # Use the larger dimension
                x_coords = np.linspace(xmin, xmax, grid_size)
                y_coords = np.linspace(ymin, ymax, grid_size)
            else:  # Custom ray trace irradiance
                D = irrad_data['D']
                w = 2 * D
                x_coords = np.linspace(-D, D, grid_size)
                y_coords = np.linspace(-D, D, grid_size)
            
            # Resize irradiance map to standard size for analysis
            try:
                from skimage.transform import resize
                z = resize(irrad_map * irrad_scal, [201, 201])
                analysis_grid_size = 201
            except ImportError:
                print("Warning: skimage not available, using original grid size")
                z = irrad_map * irrad_scal
                analysis_grid_size = grid_size
            
            # Create coordinate grids for analysis
            r = np.linspace(0, w/2, sz)
            x_analysis = np.linspace(-w/2, w/2, analysis_grid_size)
            xp, yp = np.meshgrid(x_analysis, x_analysis)
            R = np.sqrt(np.square(xp) + np.square(yp))
            
            # Calculate encircled energy
            enen = np.zeros(sz)
            for i in range(sz):
                enen[i] = np.sum(z[R <= r[i]])
            
            # Normalize encircled energy
            enen_normalized = enen / np.max(enen) if np.max(enen) > 0 else enen
            
            # Store results
            encircled_energy_data[config_key] = {
                'radius': r,
                'encircled_energy': enen_normalized,
                'total_energy': np.max(enen)
            }
            
            # Save data if requested
            if save:
                self._save_encircled_energy(config_key, r, enen_normalized)
            
            # Plot if requested
            if plot:
                self._plot_encircled_energy(config_key, r, enen_normalized)
        
        return encircled_energy_data
    
    def _save_encircled_energy(self, config_key, radius, encircled_energy):
        """Save encircled energy data to CSV"""
        data = np.column_stack((radius, encircled_energy))
        
        # Create appropriate filename
        if isinstance(config_key, str):
            filename = os.path.join(self.output_directory, 'irrad', f'encircled_energy_{config_key}.csv')
        else:
            filename = os.path.join(self.output_directory, 'irrad', f'encircled_energy_{config_key + 1}.csv')
        
        header = "Radius (mm),Encircled Energy (normalized)"
        np.savetxt(filename, data, delimiter=',', fmt='%3.6f', header=header, comments='')
        print(f"Saved encircled energy data: {filename}")
    
    def _plot_encircled_energy(self, config_key, radius, encircled_energy, dpi=300):
        """Plot and save encircled energy"""
        plt.figure(figsize=(8, 8), dpi=dpi)
        plt.rcParams['figure.dpi'] = dpi
        
        # Create title based on config key
        if isinstance(config_key, str):
            title = f"Encircled Energy - {config_key.replace('_', ' ').title()}"
        else:
            title = f"Encircled Energy - Config {config_key + 1}"
        
        plt.title(title, fontsize=16)
        plt.xlabel("Radius (mm)", fontsize=14)
        plt.ylabel("Encircled Energy", fontsize=14)
        plt.plot(radius, encircled_energy, color="black", linewidth=2)
        plt.grid(True, alpha=0.3)
        plt.ylim(0, 1)
        
        # Add some key statistics as text
        radius_50 = np.interp(0.5, encircled_energy, radius)
        radius_90 = np.interp(0.9, encircled_energy, radius)
        plt.text(0.6, 0.3, f'50% Energy: {radius_50:.1f} mm\n90% Energy: {radius_90:.1f} mm', 
                 transform=plt.gca().transAxes, fontsize=12, 
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Save plot
        if isinstance(config_key, str):
            filename = os.path.join(self.output_directory, 'irrad', f'encircled_energy_{config_key}.png')
        else:
            filename = os.path.join(self.output_directory, 'irrad', f'encircled_energy_{config_key + 1}.png')
        
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
        print(f"Saved encircled energy plot: {filename}")
        
        plt.show()
    
    def trace_and_save(self, seq_nr, nr_rays, n_traces, config_list=None, keep_in_memory=True, rotation=False, rot_angle=1):
        """
        Perform ray tracing and save results
        
        Parameters:
        -----------
        seq_nr : int
            Sequence number to trace
        nr_rays : int
            Number of rays per trace
        n_traces : int
            Number of traces to perform
        config_list : list or None
            List of configuration numbers to process (None = all configs)
        keep_in_memory : bool
            Whether to keep ray data in memory after saving (allows immediate plotting)
        rotation : bool
            Whether to apply 1 degree rotation per trace
        """
        # Check rotation warning
        if rotation and n_traces != 359 and n_traces % 359 != 0:
            print(f"Warning: {n_traces} traces will result in partial rotation. "
                  f"Consider using 359 or multiple of 359 for complete rotation.")
        
        try:
            # Store tracing parameters
            self.nr_rays = nr_rays
            self.n_traces = n_traces
            
            # Set ray distribution count
            self.core.setRayDistributionCount1(seq_nr, nr_rays)
            
            # Get surface for ray tracing
            if self.surf is None:
                surf = self.core.getSequenceImageSurface(seq_nr)
            else:
                surf = self.surf
            
            # Determine configurations to process
            num_configs = self.core.getNrConfigs()
            if num_configs == 0:
                num_configs = 1
            if config_list is None:
                configs_to_process = list(range(num_configs))
            else:
                configs_to_process = [c for c in config_list if 0 <= c < num_configs]
            
            trace_type = "rotated " if rotation else ""
            print(f"Beginning {trace_type}ray traces for {len(configs_to_process)} configurations")
            
            # Perform ray tracing
            start_time = time.time()
            if rotation:
                self._trace_configurations_with_rotation(seq_nr, n_traces, surf, configs_to_process, keep_in_memory,rot_angle)
            else:
                self._trace_configurations(seq_nr, n_traces, surf, configs_to_process, keep_in_memory)
            end_time = time.time()
            
            execution_time = end_time - start_time
            print(f"Execution time: {execution_time:.6f} seconds")
            
            # Save metadata
            metadata = {
                'seq_nr': seq_nr,
                'nr_rays': nr_rays,
                'n_traces': n_traces,
                'total_rays_traced_per_config': nr_rays * n_traces,
                'configs_processed': configs_to_process,
                'execution_time': execution_time,
                'watts_of_source': self.watts_of_source,
                'model_file': os.path.join(self.model_folder, self.model_filename),
                'surface': surf,
                'rotated': rotation
            }
            if rotation:
                metadata['rotation_angle_per_trace'] = 1.0  # degrees
            self._save_metadata(metadata)
            
        except Exception as e:
            print(f"Error in trace_and_save: {e}")
            raise

    def _trace_configurations(self, seq_nr, n_traces, surf, configs_to_process, keep_in_memory=True):
        """Private method to perform standard ray tracing - preserves NaN values"""
        # if configs_to_process[0] == 0:
        #     #Trace sequence, and save as config


        for config_num in configs_to_process:
            print(f"Starting ray trace of Configuration {config_num + 1}")
            self.core.setConfig(config_num)
            spots_list = []
            spots_dir_list = []
            
            for i in range(n_traces):
                self.core.traceRays(seq_nr, 0, 0)
                spot = self.core.getRayPos(seq_nr, 0, 0, surf)
                spot_dir = self.core.getRayDirExit(seq_nr,0,0,surf)
                npspot = np.array(spot, copy=True)
                npspot_dir = np.array(spot_dir,copy=True)
                
                # Keep all values including NaN - do not remove them
                spots_list.append(npspot)
                spots_dir_list.append(npspot_dir)
            
            # Combine all traces for this configuration
            if spots_list:
                combined_spots = np.concatenate(spots_list, axis=1)
                combined_spots_dir = np.concatenate(spots_dir_list,axis=1)
                
                # Keep in memory if requested (with NaN values preserved)
                if keep_in_memory:
                    if False:#configs_to_process[0] == 0:
                        self.ray_data[seq_nr] = combined_spots
                    else:
                        self.ray_data[config_num] = combined_spots
                        self.ray_dir_data[config_num] = combined_spots_dir
                
                # Save to CSV (with NaN values preserved)
                if False:#configs_to_process[0] == 0:
                    filename = os.path.join(self.output_directory, 'ray_trace', f'trace_{seq_nr + 1}.csv')
                else:
                    filename = os.path.join(self.output_directory, 'ray_trace', f'trace_{config_num + 1}.csv')
                    filename_dir = os.path.join(self.output_directory, 'ray_trace', f'trace_dir_{config_num + 1}.csv')
                if combined_spots.size > 0:
                    np.savetxt(filename, combined_spots.T, delimiter=',', fmt='%3.8f')
                    print(f"Saved {filename}")
                    np.savetxt(filename_dir, combined_spots_dir.T, delimiter=',', fmt='%3.8f')
                    print(f"Saved {filename_dir}")
                    
                    # Report statistics
                    total_rays = combined_spots.shape[1]
                    valid_rays = self._count_valid_rays(combined_spots)
                    print(f"  Total rays traced: {total_rays:,}")
                    print(f"  Valid rays (non-NaN): {valid_rays:,}")
                    print(f"  Transmission efficiency: {valid_rays/total_rays:.1%}")
                else:
                    print(f"Warning: No ray data for configuration {config_num + 1}")
            
            print(f"Ray trace {config_num + 1} completed")
    
    def _trace_configurations_with_rotation(self, seq_nr, n_traces, surf, configs_to_process, keep_in_memory=True,rot_angle = 1):
        """Private method to perform ray tracing with rotation - preserves NaN values"""
        for config_num in configs_to_process:
            print(f"Starting rotated ray trace of Configuration {config_num + 1}")
            self.core.setConfig(config_num)
            spots_list = []
            
            for i in range(n_traces):
                self.core.traceRays(seq_nr, 0, 0)
                spot = self.core.getRayPos(seq_nr, 0, 0, surf)
                npspot = np.array(spot, copy=True)
                
                # Apply rotation only to non-NaN values
                if npspot.size > 0:
                    rotation_angle = i * rot_angle * np.pi / 180  # Convert to radians
                    # Create mask for non-NaN values
                    valid_mask = ~np.isnan(npspot).any(axis=0)
                    if np.any(valid_mask):
                        # Apply rotation only to valid coordinates
                        npspot[:, valid_mask] = self._rot_mat(rotation_angle) @ npspot[:, valid_mask]
                
                spots_list.append(npspot)
            
            # Combine all traces for this configuration
            if spots_list:
                combined_spots = np.concatenate(spots_list, axis=1)
                
                # Keep in memory if requested (with NaN values preserved)
                if keep_in_memory:
                    self.ray_data[config_num] = combined_spots
                
                # Save to CSV (with NaN values preserved)
                filename = os.path.join(self.output_directory, 'ray_trace', f'trace_{config_num + 1}.csv')
                if combined_spots.size > 0:
                    np.savetxt(filename, combined_spots.T, delimiter=',', fmt='%3.8f')
                    print(f"Saved {filename}")
                    
                    # Report statistics
                    total_rays = combined_spots.shape[1]
                    valid_rays = self._count_valid_rays(combined_spots)
                    print(f"  Total rays traced: {total_rays:,}")
                    print(f"  Valid rays (non-NaN): {valid_rays:,}")
                    print(f"  Transmission efficiency: {valid_rays/total_rays:.1%}")
                else:
                    print(f"Warning: No ray data for configuration {config_num + 1}")
            
            print(f"Rotated ray trace {config_num + 1} completed")
    
    def load_quadoa_traces(self, config_list=None):
        """
        Load ray trace files saved manually from Quadoa (qtrace_X.csv format)
        
        Parameters:
        -----------
        config_list : list or None
            List of configuration numbers to load (None = all available qtrace files)
        """
        ray_trace_dir = os.path.join(self.output_directory, 'ray_trace')
        
        if not os.path.exists(ray_trace_dir):
            print(f"Warning: Ray trace directory not found: {ray_trace_dir}")
            return
        
        # Find qtrace files
        if config_list is None:
            qtrace_files = [f for f in os.listdir(ray_trace_dir) 
                           if f.startswith('qtrace_') and f.endswith('.csv')]
        else:
            qtrace_files = []
            for config in config_list:
                filename = f'qtrace_{config + 1}.csv'  # Convert to 1-indexed filename
                if os.path.exists(os.path.join(ray_trace_dir, filename)):
                    qtrace_files.append(filename)
                else:
                    print(f"Warning: Quadoa trace file not found: {filename}")
        
        if not qtrace_files:
            print("Warning: No qtrace files found")
            return
        
        for filename in qtrace_files:
            filepath = os.path.join(ray_trace_dir, filename)
            
            try:
                # Extract config number from filename (qtrace_X.csv -> config X-1)
                config_num = int(filename.split('_')[1].split('.')[0]) - 1
                
                # Parse Quadoa format
                ray_data = self._parse_quadoa_file(filepath)
                
                if ray_data is not None:
                    self.ray_data[config_num] = ray_data
                    
                    # Report statistics
                    total_rays = ray_data.shape[1]
                    valid_rays = self._count_valid_rays(ray_data)
                    
                    self.n_traces = 1
                    self.nr_rays = total_rays
                    
                    print(f"Loaded Quadoa trace data for configuration {config_num + 1}")
                    print(f"  Total rays: {total_rays:,}")
                    print(f"  Valid rays (non-NaN): {valid_rays:,}")
                    print(f"  Transmission efficiency: {valid_rays/total_rays:.1%}")
                else:
                    print(f"Failed to parse {filename}")
                    
            except Exception as e:
                print(f"Error loading {filename}: {e}")
    
    def _parse_quadoa_file(self, filepath):
        """
        Parse Quadoa-saved ray trace file format
        
        File format example:
        # X-Unit is m # Y-Unit is m # Field 1 / Wave 1
        3.4895343979255378e-01;2.2293844982637809e-01;
        nan;nan;
        1.0175401368369805e-01;-3.2071477232970114e-01;
        
        Parameters:
        -----------
        filepath : str
            Path to qtrace file
            
        Returns:
        --------
        ray_data : numpy.ndarray or None
            Ray positions as 2xN array (X, Y coordinates) with NaN preserved
        """
        try:
            with open(filepath, 'r') as file:
                lines = file.readlines()
            
            # Skip comment lines that start with #
            data_lines = [line.strip() for line in lines if line.strip() and not line.strip().startswith('#')]
            
            if not data_lines:
                print(f"Warning: No data found in {filepath}")
                return None
            
            # Parse coordinate pairs from semicolon-separated format
            x_coords = []
            y_coords = []
            
            for line in data_lines:
                # Split by semicolon and remove empty strings
                coords = [coord.strip() for coord in line.split(';') if coord.strip()]
                
                # Process coordinate pairs
                i = 0
                while i < len(coords):
                    if i + 1 < len(coords):
                        x_str = coords[i]
                        y_str = coords[i + 1]
                        
                        # Handle nan values
                        if x_str.lower() == 'nan':
                            x_val = np.nan
                        else:
                            try:
                                x_val = float(x_str) * 1000  # Convert from m to mm
                            except ValueError:
                                x_val = np.nan
                        
                        if y_str.lower() == 'nan':
                            y_val = np.nan
                        else:
                            try:
                                y_val = float(y_str) * 1000  # Convert from m to mm
                            except ValueError:
                                y_val = np.nan
                        
                        x_coords.append(x_val)
                        y_coords.append(y_val)
                        
                        i += 2
                    else:
                        # Odd number of coordinates - skip the last one
                        i += 1
            
            if not x_coords:
                print(f"Warning: No valid coordinate pairs found in {filepath}")
                return None
            
            # Create 2xN array
            ray_data = np.array([x_coords, y_coords])
            
            return ray_data
            
        except Exception as e:
            print(f"Error parsing Quadoa file {filepath}: {e}")
            return None
    
    def load_all_traces(self):
        """Load all ray trace CSV files from the ray_trace directory"""
        ray_trace_dir = os.path.join(self.output_directory, 'ray_trace')
        
        if not os.path.exists(ray_trace_dir):
            print(f"Warning: Ray trace directory not found: {ray_trace_dir}")
            return
        
        # Find all trace CSV files
        trace_files = [f for f in os.listdir(ray_trace_dir) 
                      if f.startswith('trace_') and f.endswith('.csv')]
        
        if not trace_files:
            print("Warning: No trace files found in ray_trace directory")
            return
        
        self._load_trace_files(ray_trace_dir, trace_files)
    
    def load_specific_traces(self, config_list):
        """
        Load specific ray trace configurations
        
        Parameters:
        -----------
        config_list : list
            List of configuration numbers to load (1-indexed for filenames)
        """
        ray_trace_dir = os.path.join(self.output_directory, 'ray_trace')
        
        if not os.path.exists(ray_trace_dir):
            print(f"Warning: Ray trace directory not found: {ray_trace_dir}")
            return
        
        # Generate expected filenames
        trace_files = []
        for config in config_list:
            filename = f'trace_{config + 1}.csv'  # Convert to 1-indexed filename
            if os.path.exists(os.path.join(ray_trace_dir, filename)):
                trace_files.append(filename)
            else:
                print(f"Warning: Trace file not found: {filename}")
        
        if trace_files:
            self._load_trace_files(ray_trace_dir, trace_files)
    
    def _load_trace_files(self, ray_trace_dir, trace_files):
        """Private method to load trace files - preserves NaN values"""
        for filename in trace_files:
            filepath = os.path.join(ray_trace_dir, filename)
            
            try:
                if os.path.getsize(filepath) > 0:
                    df = pd.read_csv(filepath, header=None)
                    ray_data = df.values.T  # Transpose to match expected format
                    
                    # Extract config number from filename (trace_X.csv -> config X-1)
                    config_num = int(filename.split('_')[1].split('.')[0]) - 1
                    self.ray_data[config_num] = ray_data
                    
                    # Report statistics
                    total_rays = ray_data.shape[1]
                    valid_rays = self._count_valid_rays(ray_data)
                    print(f"Loaded trace data for configuration {config_num + 1}")
                    print(f"  Total rays: {total_rays:,}")
                    print(f"  Valid rays (non-NaN): {valid_rays:,}")
                    print(f"  Transmission efficiency: {valid_rays/total_rays:.1%}")
                else:
                    print(f"Warning: Empty trace file: {filename}")
            except Exception as e:
                print(f"Error loading {filename}: {e}")
    
    def load_rays(self):
        """Convenience method to load all ray trace data"""
        self.load_all_traces()
    
    def plot_rays(self, figsize=(8, 8), dpi=100, alpha=0.7, color='red'):
        """
        Plot loaded ray data (removes NaN values for plotting)
        
        Parameters:
        -----------
        figsize : tuple
            Figure size (width, height)
        dpi : int
            Figure DPI
        alpha : float
            Point transparency
        color : str
            Point color
        """
        if not self.ray_data:
            print("No ray data loaded. Use load_rays() or load_quadoa_traces() first.")
            return
        
        plt.figure(figsize=figsize, dpi=dpi)
        plt.rcParams.update({'font.size': 20})
        
        for config_num, ray_positions in self.ray_data.items():
            if ray_positions.size > 0:
                # Remove NaN values for plotting
                clean_positions = self._remove_nan_values(ray_positions)
                if clean_positions.size > 0:
                    plt.scatter(clean_positions[0, :], clean_positions[1, :], 
                               s=1, color=color, alpha=alpha, 
                               label=f'Config {config_num + 1}')
        
        plt.xlabel('X (mm)')
        plt.ylabel('Y (mm)')
        plt.title('Ray Trace Results', fontsize=24)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.grid(True)
        
        # Save plot
        plot_filename = os.path.join(self.output_directory, 'ray_trace', 'ray_plot.png')
        plt.savefig(plot_filename, dpi=dpi, bbox_inches='tight')
        print(f"Ray plot saved: {plot_filename}")
        
        plt.show()
    
    def generate_irradiance_plots(self, grid_size=64, D=None, optimize_centering=False, 
                                 combine_all=False, vmin=None, vmax=None, plot_type = 'imshow', levels = None, overlay = False):
        """
        Generate irradiance plots from loaded ray data
        
        Parameters:
        -----------
        grid_size : int
            Grid resolution for irradiance calculation (default: 64)
        D : float or None
            Semi-diameter for plot bounds (mm). If None, calculated from data
        optimize_centering : bool
            Whether to perform dy optimization for annular energy balance
        combine_all : bool
            Whether to combine all configurations into a single irradiance plot
        vmin : float or None
            Minimum value for colorbar scale (default: 0)
        vmax : float or None
            Maximum value for colorbar scale (default: auto-calculated from data)
        """
        if not self.ray_data:
            print("No ray data loaded. Use load_rays() or load_quadoa_traces() first.")
            return
        
        if combine_all:
            # Combine all ray data and create single plot
            self._generate_combined_irradiance_plot(grid_size, D, optimize_centering, vmin, vmax,plot_type,levels,overlay)
        else:
            # Generate individual plots for each configuration
            self._generate_individual_irradiance_plots(grid_size, D, optimize_centering, vmin, vmax,plot_type,levels,overlay)
    
    def _generate_combined_irradiance_plot(self, grid_size, D, optimize_centering, vmin, vmax, plot_type, levels, overlay):
        """Generate a single combined irradiance plot from all configurations"""
        print("Generating combined irradiance plot from all configurations")
        
        # Check if we need ray tracing parameters (only if data came from trace_and_save)
        need_trace_params = any(config_num in self.ray_data for config_num in range(10))  # Check if we have traced data
        
        if need_trace_params and (self.nr_rays is None or self.n_traces is None):
            print("Warning: Ray tracing parameters not found. Using default scaling.")
            print("For accurate irradiance values, run trace_and_save() or provide nr_rays/n_traces when loading.")
            use_simple_scaling = True
        else:
            use_simple_scaling = False
        
        # Combine all ray data (remove NaN values for analysis)
        all_rays_x = []
        all_rays_y = []
        total_rays_traced = 0
        
        for config_num, ray_positions in self.ray_data.items():
            if ray_positions.size == 0:
                print(f"Warning: No ray data for configuration {config_num + 1}")
                continue
                
            # Count total rays (including NaN)
            total_rays_traced += ray_positions.shape[1]
            
            # Remove NaN for analysis
            clean_positions = self._remove_nan_values(ray_positions)
            if clean_positions.size > 0:
                all_rays_x.extend(clean_positions[0, :])
                all_rays_y.extend(clean_positions[1, :])
        
        if not all_rays_x:
            print("Warning: No valid ray data found")
            return
        
        x = np.array(all_rays_x)
        y = np.array(all_rays_y)
        
        # Calculate diameter if not provided
        if D is None:
            rad_dist = np.sqrt(x**2 + y**2)
            D_calc = rad_dist.max()
            print(f"Calculated diameter: {D_calc:.2f} mm for combined plot")
        else:
            D_calc = D
        
        # Generate irradiance map
        irrad_map, xbar, ybar, dy = self._generate_irradiance_map(
            x, y, grid_size, D_calc, optimize_centering
        )
        
        # Calculate irradiance scaling
        if use_simple_scaling:
            # Simple scaling without transmission efficiency correction
            num_configs = len(self.ray_data)
            valid_rays = len(x)
            input_watts = self.watts_of_source * num_configs
            bin_area_mm2 = (2 * D_calc / grid_size)**2
            bin_area_m2 = bin_area_mm2 / 1000000
            irrad_scal = input_watts / valid_rays / bin_area_m2
        else:
            # Correct scaling with transmission efficiency
            num_configs = len(self.ray_data)
            total_rays_traced_all_configs = self.nr_rays * self.n_traces * num_configs
            rays_that_made_it_total = len(x)
            
            total_input_power = self.watts_of_source * num_configs
            power_that_made_it = (rays_that_made_it_total / total_rays_traced_all_configs) * total_input_power
            
            bin_area_mm2 = (2 * D_calc / grid_size)**2
            bin_area_m2 = bin_area_mm2 / 1000000
            irrad_scal = total_input_power / bin_area_m2 / total_rays_traced_all_configs
            
            print(f"Combined plot scaling info:")
            print(f"  Total rays traced (all configs): {total_rays_traced_all_configs:,}")
            print(f"  Rays that made it to image plane: {rays_that_made_it_total:,}")
            print(f"  Transmission efficiency: {rays_that_made_it_total/total_rays_traced_all_configs:.1%}")
            print(f"  Total input power: {total_input_power:.0f} W")
            print(f"  Power that made it: {power_that_made_it:.0f} W")
            print(f"  Power per ray: {power_that_made_it/rays_that_made_it_total:.3f} W")
        
        # Store irradiance map for encircled energy analysis
        if not hasattr(self, 'irradiance_maps'):
            self.irradiance_maps = {}
        
        self.irradiance_maps["combined"] = {
            'irrad_map': irrad_map,
            'irrad_scal': irrad_scal,
            'D': D_calc,
            'grid_size': grid_size,
            'xbar': xbar,
            'ybar': ybar
        }
        
        # Save combined irradiance data
        self._save_irradiance_data("combined", irrad_map, irrad_scal, D_calc, grid_size)
        
        # Generate and save plot
        self._plot_irradiance_map("combined", irrad_map, irrad_scal, D_calc, grid_size,
                                 vmin=vmin, vmax=vmax, plot_type = plot_type, levels = levels, overlay = overlay)
    
    def _generate_individual_irradiance_plots(self, grid_size, D, optimize_centering, vmin, vmax, plot_type, levels, overlay):
        """Generate individual irradiance plots for each configuration"""
        # Check if we need ray tracing parameters (only if data came from trace_and_save)
        need_trace_params = True  # Assume we need them unless proven otherwise
        
        if self.nr_rays is None or self.n_traces is None:
            print("Warning: Ray tracing parameters not found. Using simple scaling.")
            print("For accurate irradiance values, run trace_and_save() or provide nr_rays/n_traces when loading.")
            use_simple_scaling = True
        else:
            use_simple_scaling = False
        
        for config_num, ray_positions in self.ray_data.items():
            if ray_positions.size == 0:
                print(f"Warning: No ray data for configuration {config_num + 1}")
                continue
            
            print(f"Generating irradiance plot for configuration {config_num + 1}")
            
            # Count total rays (including NaN)
            total_rays_this_config = ray_positions.shape[1]
            
            # Remove NaN values for analysis
            clean_positions = self._remove_nan_values(ray_positions)
            if clean_positions.size == 0:
                print(f"Warning: No valid rays for configuration {config_num + 1}")
                continue
            
            # Get ray positions
            x, y = clean_positions[0, :], clean_positions[1, :]
            valid_rays_this_config = len(x)
            
            # Calculate diameter if not provided
            if D is None:
                rad_dist = np.sqrt(x**2 + y**2)
                D_calc = rad_dist.max()
                print(f"Calculated diameter: {D_calc:.2f} mm for config {config_num + 1}")
            else:
                D_calc = D
            
            # Generate irradiance map
            irrad_map, xbar, ybar, dy = self._generate_irradiance_map(
                x, y, grid_size, D_calc, optimize_centering
            )
            
            # Calculate irradiance scaling
            if use_simple_scaling:
                # Simple scaling without transmission efficiency correction
                input_watts = self.watts_of_source
                bin_area_mm2 = (2 * D_calc / grid_size)**2
                bin_area_m2 = bin_area_mm2 / 1000000
                irrad_scal = input_watts / valid_rays_this_config / bin_area_m2
            else:
                # Correct scaling with transmission efficiency
                total_rays_traced_this_config = self.nr_rays * self.n_traces
                
                power_that_made_it = (valid_rays_this_config / total_rays_traced_this_config) * self.watts_of_source
                
                bin_area_mm2 = (2 * D_calc / grid_size)**2
                bin_area_m2 = bin_area_mm2 / 1000000
                irrad_scal = self.watts_of_source / bin_area_m2 / total_rays_traced_this_config
                
                print(f"Config {config_num + 1} scaling info:")
                print(f"  Rays traced: {total_rays_traced_this_config:,}")
                print(f"  Total rays in dataset: {total_rays_this_config:,}")
                print(f"  Rays that made it: {valid_rays_this_config:,}")
                print(f"  Transmission efficiency: {valid_rays_this_config/total_rays_traced_this_config:.1%}")
                print(f"  PV X: {x.max()-x.min():.0f} mm")
                print(f"  PV Y: {y.max()-y.min():.0f} mm")
                print(f"  Input power: {self.watts_of_source:.0f} W")
                print(f"  Power that made it: {power_that_made_it:.0f} W")
                print(f"  Power per ray: {power_that_made_it/valid_rays_this_config:.3f} W")
            
            # Store irradiance map for encircled energy analysis
            if not hasattr(self, 'irradiance_maps'):
                self.irradiance_maps = {}
                
            self.irradiance_maps[config_num] = {
                'irrad_map': irrad_map,
                'irrad_scal': irrad_scal,
                'D': D_calc,
                'grid_size': grid_size,
                'xbar': xbar,
                'ybar': ybar
            }
            
            # Save irradiance data
            self._save_irradiance_data(config_num, irrad_map, irrad_scal, D_calc, grid_size)
            
            # Generate and save plot
            self._plot_irradiance_map(config_num, irrad_map, irrad_scal, D_calc, grid_size,
                                     vmin=vmin, vmax=vmax, plot_type = plot_type, levels = levels, overlay = overlay)
    
    def _generate_irradiance_map(self, x, y, grid_size, D, optimize_centering=False):
        """Private method to generate irradiance map from ray positions (NaN already removed)"""
        dy = 0
        xbar = np.mean(x)
        ybar = np.mean(y)
        
        if optimize_centering:
            # Perform dy optimization (simplified version of the annular energy balance)
            radii = np.array([390, 553, 677, 782, 874]) * 250 / 874
            step = 0.25 / 4
            tolerance = 0.02
            
            while True:
                dist = np.sqrt(x**2 + (y - dy)**2)
                annular_energy = np.zeros(len(radii))
                
                for i in range(len(radii)):
                    if i == 0:
                        filt = dist < radii[i]
                    else:
                        filt = (dist < radii[i]) & (dist >= radii[i-1])
                    annular_energy[i] = np.sum(filt)
                
                first, last = annular_energy[0], annular_energy[-1]
                if first == 0:
                    step = -step
                if np.isclose(first, last, rtol=tolerance):
                    break
                dy = dy + step
                
                # Safety break to prevent infinite loop
                if abs(dy) > D:
                    print("Warning: dy optimization exceeded diameter bounds")
                    break
        
        # Apply centering
        x_centered = x# - xbar
        y_centered = (y - dy)# - ybar
        
        # Generate 2D histogram
        bin_edges_x = np.linspace(-D, D, grid_size + 1)
        bin_edges_y = np.linspace(-D, D, grid_size + 1)
        irrad_map, _, _ = np.histogram2d(y_centered, x_centered, 
                                       bins=[bin_edges_y, bin_edges_x])
        
        return irrad_map, xbar, ybar, dy
    
    def _save_irradiance_data(self, config_num, irrad_map, irrad_scal, D, grid_size):
        """Private method to save irradiance data"""
        # Create coordinate grids
        xi = np.arange(-grid_size / 2 + 0.5, grid_size / 2 + 0.5, 1)
        X, Y = np.meshgrid(xi, xi)
        X *= D / (grid_size / 2)
        Y *= D / (grid_size / 2)
        
        # Prepare data for CSV
        temp_irrad = np.round(irrad_scal * irrad_map.reshape(-1, 1), 1)
        temp_X = np.round(X.reshape(-1, 1), 1) / 1000  # Convert to meters
        temp_Y = np.round(Y.reshape(-1, 1), 1) / 1000  # Convert to meters
        temp_Z = np.zeros((len(temp_Y), 1))
        
        data = np.column_stack((temp_irrad, temp_X, temp_Y, temp_Z))
        
        # Save CSV file with appropriate filename
        if config_num == "combined":
            filename = os.path.join(self.output_directory, 'irrad', f'irrad_combined.csv')
        else:
            filename = os.path.join(self.output_directory, 'irrad', f'irrad_{config_num + 1}.csv')
        
        header = "Flux (W/m2),X (m),Y (m),Z (m)"
        np.savetxt(filename, data, delimiter=',', fmt='%3.8f', header=header, comments='')
        print(f"Saved irradiance data: {filename}")
    
    def _plot_irradiance_map(self, config_num, irrad_map, irrad_scal, D, grid_size,
                        dpi=300, vmin=None, vmax=None,
                        plot_type='imshow',   # 'imshow' (default), 'contourf', 'contour'
                        levels=None,          # None or int or sequence of levels for contours
                        overlay=False,        # if True and plot_type!='imshow', overlay contours on imshow
                        show_colorbar=True):
        """Private method to plot and save irradiance map.

        plot_type: 'imshow' (default), 'contourf', or 'contour'
        levels: None (auto), integer (number of levels), or sequence of contour levels
        overlay: if True and plot_type is 'contour' or 'contourf', draw contours over an imshow
        show_colorbar: whether to include colorbar (for imshow/contourf)
        """
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib import cm
        from matplotlib.colors import ListedColormap
        import os

        # Robust colormap sampling (works if cmap has no .colors attribute)
        base_cmap = cm.get_cmap('magma', 128)
        try:
            full_colors = base_cmap.colors
        except Exception:
            # sample colors
            full_colors = base_cmap(np.linspace(0, 1, 128))

        # build the custom colormap similar to your original intent
        list1 = full_colors[30:31]
        list2 = full_colors[45:127]
        colormap = ListedColormap(np.concatenate((list1, list2), axis=0))

        # Create figure
        plt.rc('font',size=12)
        fig, ax = plt.subplots(figsize=(8, 8), dpi=dpi)

        # Scale irradiance to kW/m^2
        scaled_irrad = irrad_scal * irrad_map / 1000.0

        # Set colorbar limits
        if vmin is None:
            vmin = 0.0
        if vmax is None:
            vmax = float(np.nanmax(scaled_irrad))

        # Prepare coordinate grid for contour plotting
        # grid_size: expected integer number of pixels / points across (matching irrad_map shape)
        ny, nx = np.array(scaled_irrad).shape
        # If user provided grid_size, use it to make symmetric grid; else infer from data shape
        if grid_size is None:
            grid_x = np.linspace(-D, D, nx)
            grid_y = np.linspace(-D, D, ny)
        else:
            grid_x = np.linspace(-D, D, grid_size)
            grid_y = np.linspace(-D, D, grid_size)

        X, Y = np.meshgrid(grid_x, grid_y)

        # Decide plotting mode
        if plot_type == 'imshow' or (overlay and plot_type in ('contour', 'contourf')):
            # Use extent so axes are in physical units
            img = ax.imshow(scaled_irrad, origin='lower', cmap=colormap,
                            vmin=vmin, vmax=vmax, extent=[-D, D, -D, D], aspect='equal')
        else:
            img = None

        # Determine contour levels
        if levels is None:
            # choose a reasonable default number of levels
            n_levels = 10
            levels = np.linspace(vmin, vmax, n_levels)
        elif isinstance(levels, int):
            levels = np.linspace(vmin, vmax, levels)
        else:
            # assume user passed a sequence
            levels = np.asarray(levels)

        # Add contours if requested
        if plot_type == 'contour':
            # contour lines only (optionally overlay on imshow)
            cs = ax.contour(X, Y, scaled_irrad, levels=levels, linewidths=0.8)
            ax.clabel(cs, inline=True, fmt='%.2f', fontsize=8)
        elif plot_type == 'contourf':
            # filled contours (optionally overlay lines on top)
            cf = ax.contourf(X, Y, scaled_irrad, levels=levels, cmap=colormap, vmin=vmin, vmax=vmax)
            if overlay:
                cs = ax.contour(X, Y, scaled_irrad, levels=levels, colors='k', linewidths=0.5)
                ax.clabel(cs, inline=True, fmt='%.2f', fontsize=7)

        # Axis labels and limits
        ax.set_xlim(-D, D)
        ax.set_ylim(-D, D)
        ax.set_xlabel('X (mm)')
        ax.set_ylabel('Y (mm)')
        ax.set_aspect('equal', adjustable='box')

        # Colorbar handling: show for imshow and contourf (or if user forces it)
        if show_colorbar and (img is not None or plot_type == 'contourf'):
            # prefer the mappable object: img (imshow) or cf (contourf)
            mappable = img if img is not None else cf
            if mappable is not None:
                cbar = fig.colorbar(mappable, ax=ax)
                cbar.set_label('Irradiance (kW/m²)', fontsize=12, labelpad=10)

        # Calculate total power (kept your original logic)
        power_total = float(np.nansum(scaled_irrad)) * self.watts_of_source / 1000.0  # result in kW

        # Save plot with appropriate filename
        suffix = plot_type
        if config_num == "combined":
            filename = os.path.join(self.output_directory, 'irrad',
                                    f'irrad_combined_{suffix}_P_{power_total:.1f}kW.png')
        else:
            filename = os.path.join(self.output_directory, 'irrad',
                                    f'irrad_{config_num + 1}_{suffix}_P_{power_total:.1f}kW.png')

        # Ensure output folder exists
        os.makedirs(os.path.dirname(filename), exist_ok=True)

        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
        print(f"Saved irradiance plot: {filename}")

        plt.show()

        # return fig, ax, filename (useful for tests / further manipulation)
        return fig, ax, filename

    def generate_1d_irrad_plot(self, integrate_dimension="x", config_key=None):
        """
        Generate 1D integrated irradiance plots from stored irradiance maps
        
        Parameters:
        -----------
        integrate_dimension : str
            Direction to integrate: "x" (sum along X to get Y profile) or "y" (sum along Y to get X profile)
        config_key : str or None
            Specific configuration to plot. If None, plots all available configurations.
            Use "combined" for combined plot, or config numbers as strings for individual configs.
        """
        # Check if irradiance maps exist
        if not hasattr(self, 'irradiance_maps') or not self.irradiance_maps:
            print("No irradiance maps found. Run generate_irradiance_plots() first.")
            return
        
        # Validate integrate_dimension
        if integrate_dimension not in ["x", "y"]:
            print("integrate_dimension must be 'x' or 'y'")
            return
        
        # Determine which configurations to plot
        if config_key is not None:
            if config_key not in self.irradiance_maps:
                print(f"Configuration '{config_key}' not found in irradiance maps.")
                print(f"Available configurations: {list(self.irradiance_maps.keys())}")
                return
            configs_to_plot = [config_key]
        else:
            configs_to_plot = list(self.irradiance_maps.keys())
        
        # Set integration parameters
        if integrate_dimension == "x":
            sum_axis = 1  # Sum along X dimension to get Y profile
            coord_label = "Y (mm)"
            file_suffix = "integrated_x"
        else:  # integrate_dimension == "y"
            sum_axis = 0  # Sum along Y dimension to get X profile
            coord_label = "X (mm)"
            file_suffix = "integrated_y"
        
        print(f"Generating 1D irradiance plots (integrating along {integrate_dimension.upper()} dimension)")
        
        # Process each configuration
        for config_key in configs_to_plot:
            irrad_data = self.irradiance_maps[config_key]
            irrad_map = irrad_data['irrad_map']
            irrad_scal = irrad_data['irrad_scal']
            D_calc = irrad_data['D']
            grid_size = irrad_data['grid_size']
            
            # Calculate integrated flux
            integrated_flux = np.sum(irrad_map, axis=sum_axis) * irrad_scal * D_calc / 1000 / grid_size
            
            # Create coordinate array from -D_calc to +D_calc
            coordinates = np.linspace(-D_calc, D_calc, grid_size)
            
            # Create plot
            plt.rc('font',size=14)
            plt.figure(figsize=(8, 8))
            plt.plot(coordinates, integrated_flux, 'b-', linewidth=2)
            plt.xlabel(coord_label)
            plt.ylabel('Integrated Flux (W/m)')
            plt.title(f'1D Integrated Irradiance - Configuration: {config_key}')
            plt.grid(True, alpha=0.3)
            
            # Save plot as PNG
            plot_filename = f"{config_key}_{file_suffix}.png"
            plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
            print(f"Saved plot: {plot_filename}")
            plt.show()
            
            # Prepare data for CSV
            csv_data = np.column_stack((coordinates, integrated_flux))
            csv_filename = f"{config_key}_{file_suffix}.csv"
            
            # Save as CSV
            header = f"{coord_label.replace(' (mm)', '')}_mm,Integrated_Flux_W_per_m"
            np.savetxt(csv_filename, csv_data, delimiter=',', header=header, comments='')
            print(f"Saved data: {csv_filename}")
            
            # Print some statistics
            max_flux = np.max(integrated_flux)
            max_position = coordinates[np.argmax(integrated_flux)]
            total_power = np.trapz(integrated_flux, coordinates) / 1000  # Convert mm to m for integration
            
            print(f"Configuration {config_key} statistics:")
            print(f"  Peak flux: {max_flux:.2e} W/m at {coord_label.split()[0]} = {max_position:.2f} mm")
            print(f"  Total integrated power: {total_power:.2f} W")
            print()