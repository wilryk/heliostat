import numpy as np

def heliostat_orientation(receiver_pos, mirror_pos, solar_az_deg, solar_el_deg, cyl_offset=0):
    """
    Compute heliostat orientation angles to reflect sunlight toward a receiver.
    """

    receiver_pos = np.array(receiver_pos, dtype=float)
    mirror_pos = np.array(mirror_pos, dtype=float)

    # Convert degrees → radians
    solar_az = np.deg2rad(solar_az_deg)
    solar_el = np.deg2rad(solar_el_deg)

    # Sun direction (from mirror toward sun)
    Rs = np.array([
        np.cos(solar_el) * np.cos(np.pi/2-solar_az),
        np.cos(solar_el) * np.sin(np.pi/2-solar_az),
        np.sin(solar_el)
    ])
    Rs /= np.linalg.norm(Rs)

    # Target direction (mirror → receiver)
    Rt = receiver_pos - mirror_pos
    focal_length = np.linalg.norm(Rt)
    Rt /= focal_length

    # Mirror normal: bisector between Rs and Rt
    n = Rs + Rt
    n /= np.linalg.norm(n)

    # Mirror normal azimuth/elevation
    rot_el = np.arcsin(n[2])
    rot_az = np.arctan2(n[1], n[0])

    # Angle of incidence
    aoi = 0.5 * np.arccos(np.clip(np.dot(Rs, Rt), -1.0, 1.0))

    # Local mirror coordinate system
    # Choose "up" vector for stability
    up = np.array([0, 0, 1])

    # Tangential (u) and sagittal (v) unit vectors on mirror plane
    u = np.cross(up, n)
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    v /= np.linalg.norm(v)

    # Project solar vector into mirror plane
    s1 = np.dot(Rs, u)
    s2 = np.dot(Rs, v)

    # Compute astigmatic rotation
    rot_astig = np.arctan2(s2, s1)

    # Convert to degrees
    rot_az_deg = np.rad2deg(rot_az)
    rot_el_deg = np.rad2deg(rot_el)
    rot_astig_deg = np.rad2deg(rot_astig)
    aoi_deg = np.rad2deg(aoi)

    # Calculate tangential/radial radii of curvature
    rad = focal_length*2
    rad_s = rad*np.cos(aoi)
    rad_t = rad/np.cos(aoi)

    return rot_az_deg, rot_el_deg, rot_astig_deg, rad_s, rad_t, aoi_deg,u,v


def heliostat_shape(rot_astig_deg, Rs, Rt):
    """
    Compute heliostat shape to focus sunlight on receiver.
    """
    rot_astig = np.deg2rad(rot_astig_deg)
    ct = 1/Rt
    cs = 1/Rs

    c0 = 1/8*(cs+ct)
    c3 = 1/4*(ct-cs)*np.sin(2*rot_astig)
    c4 = 1/8*(ct+cs)
    c5 = 1/4*(ct-cs)*np.cos(2*rot_astig)


    return c0,c3,c4,c5

def axicon_heliostat_shape_correction(u,v,sagg_vector,f_dist,f_dist_s,aoi):
    df = f_dist_s - f_dist
    dphi = 1/f_dist_s - 1/f_dist
    f_tan = 1e20
    f_sag = 1/dphi

    rad_s = f_sag*2*np.cos(aoi)
    rad_t = f_tan*2
    # Project sagittal vector into mirror plane
    s1 = np.dot(sagg_vector, u)
    s2 = np.dot(sagg_vector, v)

    # Compute astigmatic rotation
    rot_astig = np.arctan2(s2, s1)
    rot_astig_deg = np.rad2deg(rot_astig)
    # print(rot_astig*180/np.pi)

    (c0,c3,c4,c5)=heliostat_shape(rot_astig_deg,
    rad_s,
    rad_t,
    )

    return c0,c3,c4,c5


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

def get_heliostat_axicon_shape(xpos,ypos,solar_az_deg,solar_el_deg,secondary_height,receiver_height,axicon_angle_deg):
    h = secondary_height - receiver_height

    heli_field_rad = np.sqrt(xpos**2 + ypos**2)

    (
        receiver_radial_offset,
        receiver_height_offset,
        axicon_radial_intersection,
        axicon_height_intersection,
    ) = receiver_correction(
        mirror_radial_position=heli_field_rad,
        axicon_height=secondary_height,
        receiver_offset=h,
        axicon_angle_deg=axicon_angle_deg,
    )

    x_offset = xpos/heli_field_rad*receiver_radial_offset
    y_offset = ypos/heli_field_rad*receiver_radial_offset

    receiver_pos = np.array([x_offset,y_offset,secondary_height+receiver_height_offset],dtype=float)
    mirror_pos = np.array([xpos,ypos,0],dtype=float)


    dist = np.sqrt(axicon_radial_intersection**2 + axicon_height_intersection**2)
    alpha = np.deg2rad(axicon_angle_deg)
    s_prime = -np.sqrt( (h+axicon_height_intersection)**2 + axicon_radial_intersection**2)

    Rs_axicon = dist/np.tan(np.deg2rad(axicon_angle_deg))


    axicon_aoi = np.rad2deg(np.atan2(axicon_radial_intersection,h+axicon_height_intersection)+alpha)

    s = 1/(2*np.cos(np.deg2rad(axicon_aoi))/Rs_axicon - 1/s_prime)

    focal_pos = receiver_pos - mirror_pos
    focal_dist = np.sqrt(np.dot(focal_pos,focal_pos))
    focal_dist_s = focal_dist + (s + s_prime)

    (rot_az_deg,
    rot_el_deg,
    rot_astig_deg,
    rad_s,
    rad_t,
    aoi_deg,
    u_mirror,
    v_mirror,
    )=heliostat_orientation(
        receiver_pos,
        mirror_pos,
        solar_az_deg,
        solar_el_deg,
    )

    #Rotation of astigmatism for axicon is different from heliostats. How to rotate?
    (c0,c3,c4,c5)=heliostat_shape(rot_astig_deg,
    rad_s,
    rad_t,
    )

    focus_xy_u = focal_pos/focal_dist
    focus_xy_u[2] = 0
    focus_xy_u = focus_xy_u/np.linalg.norm(focus_xy_u)
    up_u = np.array((0,0,1))
    saggital_vector = np.cross(focus_xy_u,up_u)


    phi_t = 1/focal_dist
    phi_s = 1/focal_dist_s
    delta_phi_s = phi_s - phi_t

    (c0_c,c3_c,c4_c,c5_c) = axicon_heliostat_shape_correction(
        u_mirror,
        v_mirror,
        saggital_vector,
        focal_dist,
        focal_dist_s,
        np.deg2rad(aoi_deg),
    )

    qc0 = -c0
    qc3=c3/np.sqrt(6)
    qc4=-c4/np.sqrt(3)
    qc5=-c5/np.sqrt(6)

    qc3_c = -c3_c/np.sqrt(6)
    qc4_c=-c4_c/np.sqrt(3)
    qc5_c=c5_c/np.sqrt(6)

    qc3 +=qc3_c
    qc4 +=qc4_c
    qc5 +=qc5_c

    return rot_az_deg, rot_el_deg, qc3, qc4, qc5