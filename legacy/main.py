#%%
from quadoa_tools.QuadoaIrradiance import QuadoaIrradiance

ModelFolder = "C:/gitlab/heliostats"  # updated for repo relocation
filename = "heliostat_field_model_mcfg.optx"

# Initialize the class
quad = QuadoaIrradiance(
    model_folder=os.path.join(ModelFolder),
    model_filename=filename,
    output_directory=os.path.join(ModelFolder,os.path.splitext(filename)[0]),
    watts_of_source=38484.5
)



#%%
# from quadoa_tools.QuadoaIrradiance import QuadoaIrradiance
import os
from matplotlib import pyplot as plt
import numpy as np
from heliostat.noaa_solar import solar_position_calculator
from heliostat.heliostat_positions import read_helio_positions
from heliostat.heliostat_shape_solve import get_heliostat_axicon_shape
# %matplotlib qt
# load the lens file specified by user
helio_pos_name = r"heliostat\positions_selected.csv"

Year = 2026
Month = 3
Day = 21
Time = 10 / 24  # Noon, fraction of day

[solar_az_deg,solar_el_deg] = solar_position_calculator(Lat=-10.0, Lon=-52.0, TZone=None,
                                 Year=Year, Month=Month, Day=Day, Time=Time)
solar_ze_deg = 90.0 - solar_el_deg

[xpos_s,ypos_s] = read_helio_positions(os.path.join(ModelFolder,helio_pos_name))
xpos_s *=1000 #convert from m to mm
ypos_s *=1000 #convert from m to mm

secondary_height = 20000.0
receiver_height = 10000.0
axicon_angle_deg = 13.8

# quad.core.setMulticonfParam("sec_height",0,secondary_height)
# quad.core.setMulticonfParam("rec_offset",0,receiver_height - secondary_height)
# quad.core.setMulticonfParam("axi_angle",0,axicon_angle_deg)

for idx, (xpos,ypos) in enumerate(zip(xpos_s,ypos_s)):
    (rot_az_deg, rot_el_deg, qc3, qc4, qc5) = get_heliostat_axicon_shape(
        xpos,
        ypos,
        solar_az_deg,
        solar_el_deg,
        secondary_height,
        receiver_height,
        axicon_angle_deg)
    
    # quad.set_parameters(solaz=solar_az_deg,
    #                     solze=solar_ze_deg,
    #                     posx=xpos,
    #                     posy=ypos,
    #                     rot_az=rot_az_deg,
    #                     rot_el=rot_el_deg,
    #                     c3=qc3,
    #                     c4=qc4,
    #                     c5=qc5,
    #                     config_num=idx)

# filename = os.path.join(ModelFolder,"heliostat_models","config_all.optx")
# quad.core.saveModelFile(filename)
# quad.trace_and_save(3,10000,10,idx)
    


#%%
import inspect
from quadoa import *

core = QuadoaCore()

for name, member in inspect.getmembers(core, predicate=inspect.ismethod):
    if not name.startswith("_"):
        print(f"\n{name}")
        print("-" * len(name))
        print(inspect.getdoc(member))


#%%
#Sequence 4 is the receiver, 3 is the window
seq = 3
quad.trace_and_save(seq_nr= seq-1, nr_rays = 1200, n_traces = 100, rotation = False, rot_angle = 180)

quad.generate_irradiance_plots(combine_all=True,D=1000,grid_size=64,levels = 10, overlay = True,plot_type = 'contourf')

quad.calculate_encircled_energy(config_keys=["combined"])