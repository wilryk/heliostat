#%%
import os
from heliostat.noaa_solar import solar_position_calculator
from heliostat.heliostat_positions import read_helio_positions
from heliostat.heliostat_shape_solve import get_heliostat_axicon_shape

# load the heliostat position file
helio_pos_name = r"heliostat\Beam down 645 field x,y centers.csv"

#----- Fixed Geometry ------#
secondary_height = 20000.0
receiver_height = 10000.0
axicon_angle_deg = 13.8

#----- Variable time of year ------#
Year = 2026
Month = 3
Day = 21
Time = 10 / 24  # Noon, fraction of day

#----- Solar position calculation, Brazil ------#
[solar_az_deg,solar_el_deg] = solar_position_calculator(Lat=-10.0, Lon=51.9, TZone=-3,
                                 Year=Year, Month=Month, Day=Day, Time=Time)
solar_ze_deg = 90.0 - solar_el_deg

#----- Get heliostat positions in field ------#
[xpos_s,ypos_s] = read_helio_positions(os.path.join(os.getcwd(),helio_pos_name))
xpos_s *=1000 #convert from m to mm
ypos_s *=1000 #convert from m to mm



for idx, (xpos,ypos) in enumerate(zip(xpos_s,ypos_s)):
    (rot_az_deg, rot_el_deg, qc3, qc4, qc5) = get_heliostat_axicon_shape(
        xpos,
        ypos,
        solar_az_deg,
        solar_el_deg,
        secondary_height,
        receiver_height,
        axicon_angle_deg)
    print(f"Rotation Az = {rot_az_deg},Rotation El = {rot_el_deg}")
