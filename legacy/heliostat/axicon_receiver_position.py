
from numpy import pi, sin, cos, tan

def receiver_correction(mirror_radial_position,axicon_height,receiver_offset, axicon_angle_deg):
    x_r = -receiver_offset * sin(2*axicon_angle_deg*pi/180)
    y_r = receiver_offset * cos(2*axicon_angle_deg*pi/180)

    x_m = mirror_radial_position
    y_m = -axicon_height

    m = (y_r-y_m)/(x_r-x_m)

    x_a = (m*x_m-y_m)/(m-tan(axicon_angle_deg*pi/180))
    y_a = tan(axicon_angle_deg*pi/180)*x_a

    receiver_radial_offset = x_r
    receiver_height_offset = y_r

    axicon_radial_intersection = x_a
    axicon_height_intersection = y_a

    return receiver_radial_offset, receiver_height_offset, axicon_radial_intersection, axicon_height_intersection