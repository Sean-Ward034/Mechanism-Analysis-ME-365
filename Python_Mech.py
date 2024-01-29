import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# TakeHome Project Spring 2023
a_pin = 1.00
b_pin = 1.50
c_pin = 2.50
d_pin = 2.315
a_slider = 4.00
b_slider = 3.00
c_slider = 1.542635

angle_G5CE = 21.2505
g5C = 1.6094

omega2 = -4
alpha2 = 0

mass_per_unit_length = 3  #kg
m2 = a_pin*mass_per_unit_length
m3 = b_pin*mass_per_unit_length
m4 = a_slider*mass_per_unit_length
m5 = 6
m6 = 2
iG2 = (1/12)*m2*(a_pin**2)
iG3 = (1/12)*m3*(b_pin**2)
iG4 = (1/12)*m4*(a_slider**2)
iG5 = 6.0208

angle_O2O4 = 48.33
angle_slider = 30
k1 = d_pin/a_pin
k2 = d_pin/c_pin
k3 = ((a_pin**2)-(b_pin**2)+(c_pin**2)+(d_pin**2))/(2*a_pin*c_pin)
k4 = d_pin/b_pin
k5 = ((c_pin**2)-(d_pin**2)-(a_pin**2)-(b_pin**2))/(2*a_pin*b_pin)

d_min = -5.93009
d_max = -1.68001
theta2_for_dmin = 166
theta2_for_dmax = 302
F = 5000
g = 9.81

def positionanalysis_PinJointed (theta2_pin):
    theta3_and_theta4 = ["", ""]
    theta2_radians = math.radians(theta2_pin)
    aa = math.cos(theta2_radians) - k1 - k2 * math.cos(theta2_radians) + k3
    bb = -2*math.sin(theta2_radians)
    cc = k1-(k2+1)*math.cos(theta2_radians)+k3
    dd = math.cos(theta2_radians) - k1 + (k4 * math.cos(theta2_radians)) + k5
    ee = bb
    ff = k1+(k4-1)*math.cos(theta2_radians)+k5
    theta3_rad = 2 * math.atan((-ee + ((ee ** 2 - 4 * dd * ff)**0.5)) / (2 * dd))
    theta3_deg = math.degrees(theta3_rad)

    theta4_rad = 2 * math.atan((-bb + (bb ** 2 - 4 * aa * cc)**0.5) / (2 * aa))
    theta4_deg = math.degrees(theta4_rad)

    theta3_and_theta4 [0] = theta3_deg
    theta3_and_theta4[1] = theta4_deg
    return theta3_and_theta4

def positionanalysis_SliderCrank (theta2_slider):
    theta3_and_d = ["", ""]
    theta2_slider_rad = math.radians(theta2_slider)
    theta3_slider_rad = math.asin((a_slider*math.sin(theta2_slider_rad)-c_slider)/b_slider)
    theta3_slider_deg = math.degrees(theta3_slider_rad)
    d = a_slider*math.cos(theta2_slider_rad)-b_slider*math.cos(theta3_slider_rad)
#    print (f'theta3_slider_deg = {theta3_slider_deg}, d = {d}')
    theta3_and_d [0] = theta3_slider_deg
    theta3_and_d[1] = d
    return theta3_and_d


def velocityAnalysis_pinJointed (theta2, theta3, theta4, omega2):
    omega3_and_omega4 = ["", ""]
    angle1_rad = math.radians(theta4 - theta2)
    angle2_rad = math.radians(theta3 - theta4)
    angle3_rad = math.radians(theta2 - theta3)
    angle4_rad = math.radians(theta4 -theta3)
    omega3 = a_pin * omega2 * math.sin(angle1_rad) / (b_pin * math.sin(angle2_rad))
    omega4 = (a_pin*omega2*math.sin(angle3_rad))/(c_pin*math.sin(angle4_rad))
    omega3_and_omega4 [0] = omega3
    omega3_and_omega4[1] = omega4
    return omega3_and_omega4

def velocityAnalysis_sliderCrank (theta4, theta5, omega4):
    omega5_and_d_dot = ["", ""]
    omega5 = a_slider*math.cos(math.radians(theta4))*omega4/(b_slider*math.cos(math.radians(theta5)))
    d_dot = -a_slider*omega4*math.sin(math.radians(theta4))+b_slider*omega5*math.sin(math.radians(theta5))
    omega5_and_d_dot [0] = omega5
    omega5_and_d_dot[1] = d_dot
    return omega5_and_d_dot

def accelAnalysis_pinJointed(theta2, theta3, theta4, omega2, omega3, omega4, alpha2):
    alpha3_and_alpha4 = ["", ""]
    theta2_rad = math.radians(theta2)
    theta3_rad = math.radians(theta3)
    theta4_rad = math.radians(theta4)
    aA = c_pin*math.sin(theta4_rad)
    bB = b_pin*math.sin(theta3_rad)
    cC1 = a_pin*alpha2*math.sin(theta2_rad) + a_pin*(omega2**2)*math.cos(theta2_rad)
    cC2 = b_pin*(omega3**2)*math.cos(theta3_rad) - c_pin*(omega4**2)*math.cos(theta4_rad)
    cC = cC1 + cC2
    dD = c_pin*math.cos(theta4_rad)
    eE = b_pin*math.cos(theta3_rad)
    fF1 = a_pin*alpha2*math.cos(theta2_rad) - a_pin*(omega2**2)*math.sin(theta2_rad)
    fF2 = -b_pin*(omega3**2)*math.sin(theta3_rad) + c_pin*(omega4**2)*math.sin(theta4_rad)
    fF = fF1 + fF2
    alpha3 = (cC*dD - aA*fF)/(aA*eE - bB*dD)
    alpha4 = (cC*eE - bB*fF)/(aA*eE - bB*dD)
    alpha3_and_alpha4 [0] = alpha3
    alpha3_and_alpha4[1] = alpha4
    return alpha3_and_alpha4

def accelAnalysis_sliderCrank(theta4_slider, theta5_slider, omega4, omega5):
    alpha5_and_d_dotdot = ["", ""]
    theta4_slider_rad = math.radians(theta4_slider)
    theta5_slider_rad = math.radians(theta5_slider)
    numerator1 = a_slider*alpha4*math.cos(theta4_slider_rad) - a_slider*(omega4**2)*math.sin(theta4_slider_rad)
    numerator2 = b_slider*(omega5**2)*math.sin(theta5_slider_rad)
    alpha5 = (numerator1 + numerator2)/(b_slider*math.cos(theta5_slider_rad))
    d_dotdot1 = -a_slider*alpha4*math.sin(theta4_slider_rad) - a_slider*(omega4**2)*math.cos(theta4_slider_rad)
    d_dotdot2 = b_slider*alpha5*math.sin(theta5_slider_rad) + b_slider*(omega5**2)*math.cos(theta5_slider_rad)
    d_dotdot = d_dotdot1+d_dotdot2
    alpha5_and_d_dotdot[0] = alpha5
    alpha5_and_d_dotdot[1] = d_dotdot
    return alpha5_and_d_dotdot

def accel_of_center_of_masses_pinJointed (theta2x):
    aG2x_aG2y_aG3x_aG3y_aG4x_aG4y = ['', "", "", "", "", ""]
    theta_G2O2_rad = math.radians(theta2x)
    Theta_G3A_rad = math.radians(theta3x)
    theta_G4O4_rad = math.radians(theta4x)
    g2O2 = a_pin/2
    g3A = b_pin/2
    g4O4 = a_slider/2
    a_G2x = g2O2*(-(omega2**2)*math.cos(theta_G2O2_rad) - alpha2*math.sin(theta_G2O2_rad))
    a_G2y = g2O2 * (-(omega2 ** 2) * math.sin(theta_G2O2_rad) + alpha2 * math.cos(theta_G2O2_rad))
    a_G3x = a_pin*(-(omega2**2)*math.cos(theta_G2O2_rad) - alpha2*math.sin(theta_G2O2_rad))\
            + g3A*(-(omega3**2)*math.cos(Theta_G3A_rad) - alpha3*math.sin(Theta_G3A_rad))
    a_G3y = a_pin*(-(omega2**2)*math.sin(theta_G2O2_rad) + alpha2*math.cos(theta_G2O2_rad))\
            + g3A*(-(omega3**2)*math.sin(Theta_G3A_rad) + alpha3*math.cos(Theta_G3A_rad))
    a_G4x = g4O4*(-(omega4**2)*math.cos(theta_G4O4_rad) - alpha4*math.sin(theta_G4O4_rad))
    a_G4y = g4O4 * (-(omega4 ** 2) * math.sin(theta_G4O4_rad) + alpha4 * math.cos(theta_G4O4_rad))
    aG2x_aG2y_aG3x_aG3y_aG4x_aG4y[0] = a_G2x
    aG2x_aG2y_aG3x_aG3y_aG4x_aG4y[1] = a_G2y
    aG2x_aG2y_aG3x_aG3y_aG4x_aG4y[2] = a_G3x
    aG2x_aG2y_aG3x_aG3y_aG4x_aG4y[3] = a_G3y
    aG2x_aG2y_aG3x_aG3y_aG4x_aG4y[4] = a_G4x
    aG2x_aG2y_aG3x_aG3y_aG4x_aG4y[5] = a_G4y
    return [float(a_G2x), float(a_G2y), float(a_G3x), float(a_G3y), float(a_G4x), float(a_G4y)]

def accel_of_center_of_masses_sliderCrank ():
    aG5x_aG5y = ["", ""]
    theta_CO4_rad = math.radians(theta4x)
    theta_G5C_rad = math.radians(theta5x - angle_G5CE)

    a_G5x = a_slider*((-omega4**2)*math.cos(theta_CO4_rad) - alpha4*math.sin(theta_CO4_rad))\
            +g5C*((-omega5**2)*math.cos(theta_G5C_rad) - alpha5*math.sin(theta_G5C_rad))
    a_G5y = a_slider*((-omega4**2)*math.sin(theta_CO4_rad) + alpha4*math.cos(theta_CO4_rad))\
            +g5C*((-omega5**2)*math.sin(theta_G5C_rad) + alpha5*math.cos(theta_G5C_rad))
    aG5x_aG5y[0] = a_G5x
    aG5x_aG5y[1] = a_G5y
    return [float(a_G5x), float(a_G5y)]
    
def calculate_position_vectors(theta2x, theta3x, theta4x, a_pin, b_pin, c_pin, d_pin, theta2, theta3, theta4, omega2, omega3, omega4, alpha2):
    c = 4.0
    BC = 1.5
    alpha3_and_alpha4 = ["", ""]
    theta2_rad = math.radians(theta2)
    theta3_rad = math.radians(theta3)
    theta4_rad = math.radians(theta4)
    aA = c_pin*math.sin(theta4_rad)
    bB = b_pin*math.sin(theta3_rad)
    cC1 = a_pin*alpha2*math.sin(theta2_rad) + a_pin*(omega2**2)*math.cos(theta2_rad)
    cC2 = b_pin*(omega3**2)*math.cos(theta3_rad) - c_pin*(omega4**2)*math.cos(theta4_rad)
    cC = cC1 + cC2
    # Define the angles for each R vector based on theta2x
    theta_65 = 21.25
    theta_R12 = math.radians(theta2x + 180)
    theta_R32 = math.radians(theta2x)
    theta_R23 = math.radians(theta3x - 180)
    theta_R34 = math.radians(theta4x)
    theta_R43 = math.radians(theta3x)
    theta_R14 = math.radians(theta4x + 180)
    theta_R54 = math.radians(theta4x)
    theta_R45 = math.radians(theta5x - theta_65 + 180)
    theta_R65 = math.radians(theta5x + theta_65)

    # Calculate the position vectors for R12
    R12x = (a_pin / 2) * math.cos(theta_R12)
    R12y = (a_pin / 2) * math.sin(theta_R12)

    # Calculate the position vectors for R32
    R32x = (a_pin / 2) * math.cos(theta_R32)
    R32y = (a_pin / 2) * math.sin(theta_R32)

    # Calculate the position vectors for R23
    R23x = (b_pin / 2) * math.cos(theta_R23)
    R23y = (b_pin / 2) * math.sin(theta_R23)

    # Calculate the position vectors for R43
    R43x = (b_pin / 2) * math.cos(theta_R43)
    R43y = (b_pin / 2) * math.sin(theta_R43)
    
    # Calculate the postion vectors for R34
    R34x = (((c / 2) - (BC)) * math.cos(theta_R34))
    R34y = (((c / 2) - (BC)) * math.sin(theta_R34))

    # Calculate the position vectors for R14
    R14x = (c / 2) * math.cos(theta_R14)
    R14y = (c / 2) * math.sin(theta_R14)

    # Calculate the position vectors for R54
    R54x = (c / 2) * math.cos(theta_R54)
    R54y = (c / 2) * math.sin(theta_R54)

    # Calculate the position vectors for R45
    R45x = (g5C) * math.cos(theta_R45)
    R45y = (g5C) * math.sin(theta_R45)

    # Calculate the position vectors for R65
    R65x = (g5C) * math.cos(theta_R65)
    R65y = (g5C) * math.sin(theta_R65)
    
    # Return all calculated vectors
    return {
        'R12x': R12x, 'R12y': R12y,
        'R32x': R32x, 'R32y': R32y,
        'R23x': R23x, 'R23y': R23y,
        'R43x': R43x, 'R43y': R43y,
        'R34x': R34x, 'R34y': R34y,
        'R14x': R14x, 'R14y': R14y,
        'R54x': R54x, 'R54y': R54y,
        'R45x': R45x, 'R45y': R45y,
        'R65x': R65x, 'R65y': R65y,
    }
    
def solve_forceMatrix(R12x, R12y, R32x, R32y, R23x, R23y, R43x, R43y, R34x, R34y, R14x, R14y, R54x, R54y, R45x, R45y, R65x, R65y, m2, m3, m4, m5, m6, iG2, iG3, iG4, iG5, F, g, a_G2x, a_G2y, a_G3x, a_G3y, a_G4x, a_G4y, a_G5x, a_G5y, a_G6x, a_G6y, alpha2, alpha3, alpha4, alpha5):
    force_matrix = np.array([
        [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-R12y, R12y, -R32y, R32x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, R23y, -R23x, -R43y, R43x, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, -1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, -1, 0, 1, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, R34y, -R34x, -R14y, R14x, -R54y, R54x, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, R45y, -R45x, -R65y, R65x, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, np.cos(np.radians(30)), -np.sin(np.radians(30)), 0]
    ])
    
    solutions_matrix = np.array([
       	[m2*a_G2x],
        [m2*a_G2y + m2*g],
        [iG2*alpha2],
        [m3*a_G3x],
        [m3*a_G3y + m3*g],
        [iG3*alpha3],
        [m4*a_G4x],
        [m4*a_G4y + m4*g],
        [iG4*alpha4],
        [m5*a_G5x],
        [m5*a_G5y + m5*g + F],
        [iG5*alpha5],
        [m6*a_G6x],
        [m6*a_G6y + m6*g],
        [0]
    ])
    # Solve the force matrix to get the force components
    solutions = np.linalg.solve(force_matrix, solutions_matrix)
    # Return the solutions with 'F' values in the specified format
    return {
        'F12x': solutions[0][0], 'F12y': solutions[1][0],
        'F32x': solutions[2][0], 'F32y': solutions[3][0],
        'F43x': solutions[4][0], 'F43y': solutions[5][0],
        'F14x': solutions[6][0], 'F14y': solutions[7][0],
        'F54x': solutions[8][0], 'F54y': solutions[9][0],
        'F65x': solutions[10][0], 'F65y': solutions[11][0],
        'F16x': solutions[12][0], 'F16y': solutions[13][0],
        'T12': solutions[14][0]
    }

def calculate_magnitudes(F12x, F12y, F32x, F32y, F43x, F43y, F14x, F14y, F54x, F54y, F65x, F65y, F16x, F16y, T12):
    # Calculate the magnitudes of the forces based on their x and y components
    F02 = np.sqrt(F12x**2 + F12y**2)
    FA = np.sqrt(F32x**2 + F32y**2)
    FB = np.sqrt(F43x**2 + F43y**2)
    F04 = np.sqrt(F14x**2 + F14y**2)
    FC = np.sqrt(F54x**2 + F54y**2)
    FD = np.sqrt(F65x**2 + F65y**2)
    F16 = np.sqrt(F16x**2 + F16y**2)

    # Return the magnitudes in a new dictionary
    return {
        'F02': F02,
        'FA': FA,
        'FB': FB,
        'F04': F04,
        'FC': FC,
        'FD': FD,
        'F16': F16,
        'T12': T12
    }

# Data lists for plotting
theta2x_list = []
theta3x_list = []
theta4x_list = []
theta5x_list = []
d_list = []

omega2_list = []
omega3_list = []
omega4_list = []
omega5_list = []
d_dot_list = []

alpha2_list = []
alpha3_list = []
alpha4_list = []
alpha5_list = []
d_dotdot_list = []

aG2x_list = []
aG2y_list = []
aG3x_list = []
aG3y_list = []
aG4x_list = []
aG4y_list = []
aG5x_list = []
aG5y_list = []
aG6x_list = []
aG6y_list = []

R12x_list = []
R12y_list = []
R32x_list = []
R32y_list = []
R23x_list = []
R23y_list = []
R43x_list = []
R43y_list = []
R34x_list = []
R34y_list = []
R14x_list = []
R14y_list = []
R54x_list = []
R54y_list = []
R45x_list = []
R45y_list = []
R65x_list = []
R65y_list = []

FO2_list = []
FA_list = []
FB_list = []
FO4_list = []
FC_list = []
FD_list = []
T12_list = []

# Initialize variables to record the maximum magnitudes and corresponding angles
max_F02, angle_max_F02 = 0, 0
max_FA, angle_max_FA = 0, 0
max_FB, angle_max_FB = 0, 0
max_F04, angle_max_F04 = 0, 0
max_FC, angle_max_FC = 0, 0
max_FD, angle_max_FD = 0, 0
max_T12, angle_max_T12 = 0, 0

# MAIN LOOP to calculate the parameters for each angle theta2x from 0 to 360 degrees

for x in range(0, 361, 1):
    theta2x = x
    theta2 = theta2x + 180 - angle_O2O4
    theta3_and_theta4 = positionanalysis_PinJointed(theta2)
    theta3 = theta3_and_theta4[0]
    theta3x = theta3 - 180 + angle_O2O4
    if theta3x < 0:
        theta3x = theta3x + 360

    theta4 = theta3_and_theta4[1]
    theta4x = theta4 - 180 + angle_O2O4
    if theta4x < 0:
        theta4x = theta4x + 360

    gamma = abs(theta3x - theta4x)

    theta4_slider = theta4 - 180 + angle_O2O4 + angle_slider
    theta5_and_d = positionanalysis_SliderCrank(theta4_slider)
    theta5_slider = theta5_and_d[0]
    theta5x = 180 - angle_slider + theta5_slider
    if theta5x < 0:
        theta5x = theta5x + 360

    d = theta5_and_d[1]

    theta2x_list.append(theta2x)
    theta3x_list.append(theta3x)
    theta4x_list.append(theta4x)
    theta5x_list.append(theta5x)
    d_list.append(d)

    omega3_and_omega4 = velocityAnalysis_pinJointed(theta2, theta3, theta4, omega2)
    omega3 = omega3_and_omega4[0]
    omega4 = omega3_and_omega4[1]
    omega5_and_d_dot = velocityAnalysis_sliderCrank (theta4_slider, theta5_slider, omega4)
    omega5 = omega5_and_d_dot[0]
    d_dot = omega5_and_d_dot[1]

    omega2_list.append(omega2)
    omega3_list.append(omega3)
    omega4_list.append(omega4)
    omega5_list.append(omega5)
    d_dot_list.append(d_dot)

#   alpha2 = calculate_alpha2(omega2, delta)
    alpha3_and_alpha4 = accelAnalysis_pinJointed(theta2, theta3, theta4, omega2, omega3, omega4, alpha2)
    alpha3 = alpha3_and_alpha4[0]
    alpha4 = alpha3_and_alpha4[1]
    alpha5_and_d_dotdot = accelAnalysis_sliderCrank(theta4_slider, theta5_slider, omega4, omega5)
    alpha5 = alpha5_and_d_dotdot[0]
    d_dotdot = alpha5_and_d_dotdot[1]
    
    alpha2_list.append(alpha2)
    alpha3_list.append(alpha3)
    alpha4_list.append(alpha4)
    alpha5_list.append(alpha5)
    d_dotdot_list.append(d_dotdot)

    aG2x_aG2y_aG3x_aG3y_aG4x_aG4y = accel_of_center_of_masses_pinJointed(theta2x)
    a_G2x = aG2x_aG2y_aG3x_aG3y_aG4x_aG4y[0]
    a_G2y = aG2x_aG2y_aG3x_aG3y_aG4x_aG4y[1]
    a_G3x = aG2x_aG2y_aG3x_aG3y_aG4x_aG4y[2]
    a_G3y = aG2x_aG2y_aG3x_aG3y_aG4x_aG4y[3]
    a_G4x = aG2x_aG2y_aG3x_aG3y_aG4x_aG4y[4]
    a_G4y = aG2x_aG2y_aG3x_aG3y_aG4x_aG4y[5]

    aG2x_list.append(a_G2x)
    aG2y_list.append(a_G2y)
    aG3x_list.append(a_G3x)
    aG3y_list.append(a_G3y)
    aG4x_list.append(a_G4x)
    aG4y_list.append(a_G4y)

    aG5x_aG5y = accel_of_center_of_masses_sliderCrank()
    a_G5x = aG5x_aG5y[0]
    a_G5y = aG5x_aG5y[1]
    aG5x_list.append(a_G5x)
    aG5y_list.append(a_G5y)

    # Calculate a_G6x and a_G6y using d_dotdot, angle_slider, and add to list
    a_G6x = d_dotdot * math.cos(math.radians(angle_slider))
    a_G6y = -d_dotdot * math.sin(math.radians(angle_slider))
    aG6x_list.append(a_G6x)
    aG6y_list.append(a_G6y)
    
    # Call the function to calculate position vectors and append to lists
    R_vectors = calculate_position_vectors(theta2x, theta3x, theta4x, a_pin, b_pin, c_pin, d_pin, theta2, theta3, theta4, omega2, omega3, omega4, alpha2)
    R12x_list.append(R_vectors['R12x'])
    R12y_list.append(R_vectors['R12y'])
    R32x_list.append(R_vectors['R32x'])
    R32y_list.append(R_vectors['R32y'])
    R23x_list.append(R_vectors['R23x'])
    R23y_list.append(R_vectors['R23y'])
    R43x_list.append(R_vectors['R43x'])
    R43y_list.append(R_vectors['R43y'])
    R34x_list.append(R_vectors['R34x'])
    R34y_list.append(R_vectors['R34y'])
    R14x_list.append(R_vectors['R14x'])
    R14y_list.append(R_vectors['R14y'])
    R54x_list.append(R_vectors['R54x'])
    R54y_list.append(R_vectors['R54y'])
    R45x_list.append(R_vectors['R45x'])
    R45y_list.append(R_vectors['R45y'])
    R65x_list.append(R_vectors['R65x'])
    R65y_list.append(R_vectors['R65y'])
    
    # Now extract each component from the returned dictionary
    R12xx = R_vectors['R12x']
    R12yy = R_vectors['R12y']
    R32xx = R_vectors['R32x']
    R32yy = R_vectors['R32y']
    R23xx = R_vectors['R23x']
    R23yy = R_vectors['R23y']
    R43xx = R_vectors['R43x']
    R43yy = R_vectors['R43y']
    R34xx = R_vectors['R34x']
    R34yy = R_vectors['R34y']
    R14xx = R_vectors['R14x']
    R14yy = R_vectors['R14y']
    R54xx = R_vectors['R54x']
    R54yy = R_vectors['R54y']
    R45xx = R_vectors['R45x']
    R45yy = R_vectors['R45y']
    R65xx = R_vectors['R65x']
    R65yy = R_vectors['R65y']
    
    # Here, call your solve_forceMatrix function with the required parameters
    force_components = solve_forceMatrix(
        R12xx, R12yy, R32xx, R32yy, R23xx, R23yy, R43xx, R43yy, R34xx, R34yy, 
        R14xx, R14yy, R54xx, R54yy, R45xx, R45yy, R65xx, R65yy, 
        m2, m3, m4, m5, m6, iG2, iG3, iG4, iG5, 
        F, g, a_G2x, a_G2y, a_G3x, a_G3y, a_G4x, a_G4y, 
        a_G5x, a_G5y, a_G6x, a_G6y, alpha2, alpha3, alpha4, alpha5
    )

    # Now, calculate the magnitudes of these force components
    magnitudes = calculate_magnitudes(
        force_components['F12x'], force_components['F12y'], force_components['F32x'], force_components['F32y'],
        force_components['F43x'], force_components['F43y'], force_components['F14x'], force_components['F14y'],
        force_components['F54x'], force_components['F54y'], force_components['F65x'], force_components['F65y'],
        force_components['F16x'], force_components['F16y'], force_components['T12']
    )

    # Append the calculated magnitudes to your lists
    FO2_list.append(magnitudes['F02'])
    FA_list.append(magnitudes['FA'])
    FB_list.append(magnitudes['FB'])
    FO4_list.append(magnitudes['F04'])
    FC_list.append(magnitudes['FC'])
    FD_list.append(magnitudes['FD'])
    T12_list.append(magnitudes['T12'])
    
    # Update the max magnitude values if current magnitudes are greater
    max_F02 = max(max_F02, magnitudes['F02'])
    max_FA = max(max_FA, magnitudes['FA'])
    max_FB = max(max_FB, magnitudes['FB'])
    max_F04 = max(max_F04, magnitudes['F04'])
    max_FC = max(max_FC, magnitudes['FC'])
    max_FD = max(max_FD, magnitudes['FD'])
    max_T12 = max(max_T12, magnitudes['T12'])
    
    # Update the max magnitude values and their corresponding angles
    if magnitudes['F02'] > max_F02:
        max_F02, angle_max_F02 = magnitudes['F02'], theta2x
    if magnitudes['FA'] > max_FA:
        max_FA, angle_max_FA = magnitudes['FA'], theta2x
    if magnitudes['FB'] > max_FB:
        max_FB, angle_max_FB = magnitudes['FB'], theta2x
    if magnitudes['F04'] > max_F04:
        max_F04, angle_max_F04 = magnitudes['F04'], theta2x
    if magnitudes['FC'] > max_FC:
        max_FC, angle_max_FC = magnitudes['FC'], theta2x
    if magnitudes['FD'] > max_FD:
        max_FD, angle_max_FD = magnitudes['FD'], theta2x
    if magnitudes['T12'] > max_T12:
        max_T12, angle_max_T12 = magnitudes['T12'], theta2x
    
# Position DataFrame
position_data = {
    'Theta2x': theta2x_list,
    'Theta3x': theta3x_list,
    'Theta4x': theta4x_list,
    'Theta5x': theta5x_list,
    'D': d_list
}
position_df = pd.DataFrame(position_data)

# Velocity DataFrame
velocity_data = {
    'Theta2x': theta2x_list,
    'Omega2': omega2_list,
    'Omega3': omega3_list,
    'Omega4': omega4_list,
    'Omega5': omega5_list,
    'D_dot': d_dot_list
}
velocity_df = pd.DataFrame(velocity_data)

# Acceleration DataFrame
acceleration_data = {
    'Theta2x': theta2x_list,
    'Alpha2': alpha2_list,
    'Alpha3': alpha3_list,
    'Alpha4': alpha4_list,
    'Alpha5': alpha5_list,
    'D_dotdot': d_dotdot_list
}
acceleration_df = pd.DataFrame(acceleration_data)

# Acceleration of Center of Masses DataFrame
acceleration_cm_data = {
    'Theta2x': theta2x_list,
    'aG2x': aG2x_list,
    'aG2y': aG2y_list,
    'aG3x': aG3x_list,
    'aG3y': aG3y_list,
    'aG4x': aG4x_list,
    'aG4y': aG4y_list,
    'aG5x': aG5x_list,
    'aG5y': aG5y_list,
    'aG6x': aG6x_list,
    'aG6y': aG6y_list
}
acceleration_cm_df = pd.DataFrame(acceleration_cm_data)

# DataFrame for the first set of position vectors
position_vectors_data_1 = {
    'Theta2x': theta2x_list,
    'R12x': R12x_list,
    'R12y': R12y_list,
    'R32x': R32x_list,
    'R32y': R32y_list,
    'R23x': R23x_list,
    'R23y': R23y_list,
    'R43x': R43x_list,  # Assuming you meant R43x instead of R34x for symmetry
    'R43y': R43y_list   # Assuming you meant R43y instead of R34y for symmetry
}
position_vectors_df_1 = pd.DataFrame(position_vectors_data_1)

# DataFrame for the second set of position vectors
position_vectors_data_2 = {
    'Theta2x': theta2x_list,
    'R14x': R14x_list,  # Assuming R14x/R14y are calculated and lists are initialized
    'R14y': R14y_list,
    'R34x': R34x_list,
    'R34y': R34y_list,
    'R54x': R54x_list,  # Assuming R54x/R54y are calculated and lists are initialized
    'R54y': R54y_list,
    'R45x': R45x_list,
    'R45y': R45y_list,
    'R65x': R65x_list,
    'R65y': R65y_list
}
position_vectors_df_2 = pd.DataFrame(position_vectors_data_2)

# Create a DataFrame for the maximum magnitudes
max_magnitudes = pd.DataFrame({
    'Max Magnitude': ['F02', 'FA', 'FB', 'F04', 'FC', 'FD', 'T12'],
    'Value': [max_F02, max_FA, max_FB, max_F04, max_FC, max_FD, max_T12],
    'Angle at Max': [angle_max_F02, angle_max_FA, angle_max_FB, angle_max_F04, angle_max_FC, angle_max_FD, angle_max_T12]
})
print("Force Analysis Data for Links")
print(max_magnitudes)

# Now display only the specific rows for the given theta2x values
print("Position Analysis Data for Theta2x = 120:")
print(position_df[position_df['Theta2x'] == 120])

print("\nVelocity Analysis Data for Theta2x = 120:")
print(velocity_df[velocity_df['Theta2x'] == 120])

print("\nAcceleration Analysis Data for Theta2x = 120:")
print(acceleration_df[acceleration_df['Theta2x'] == 120])

print("\nAcceleration of Center of Masses Data for Theta2x = 120:")
print(acceleration_cm_df[acceleration_cm_df['Theta2x'] == 120])

print("\nPosition Vectors (table 1) for Theta2x = 120:")
print(position_vectors_df_1[position_vectors_df_1['Theta2x'] == 120])

print("\nPosition Vectors (table 2) for Theta2x = 120:")
print(position_vectors_df_2[position_vectors_df_2['Theta2x'] == 120])

# Plotting the data using matplotlib
plt.figure(figsize=(15, 10))

# Plot for Position Analysis
plt.subplot(3, 2, 1)
plt.plot(position_df['Theta2x'], position_df['Theta3x'], label='Theta3x')
plt.plot(position_df['Theta2x'], position_df['Theta4x'], label='Theta4x')
plt.plot(position_df['Theta2x'], position_df['Theta5x'], label='Theta5x')
plt.xlabel('Theta2x (degrees)')
plt.ylabel('Angles (degrees)')
plt.title('Position Analysis')
plt.legend()

# Plot for Velocity Analysis
plt.subplot(3, 2, 2)
plt.plot(velocity_df['Theta2x'], velocity_df['Omega2'], label='Omega2')
plt.plot(velocity_df['Theta2x'], velocity_df['Omega3'], label='Omega3')
plt.plot(velocity_df['Theta2x'], velocity_df['Omega4'], label='Omega4')
plt.plot(velocity_df['Theta2x'], velocity_df['Omega5'], label='Omega5')
plt.xlabel('Theta2x (degrees)')
plt.ylabel('Angular Velocities (rad/s)')
plt.title('Velocity Analysis')
plt.legend()

# Plot for Acceleration Analysis
plt.subplot(3, 2, 3)
plt.plot(acceleration_df['Theta2x'], acceleration_df['Alpha2'], label='Alpha2')
plt.plot(acceleration_df['Theta2x'], acceleration_df['Alpha3'], label='Alpha3')
plt.plot(acceleration_df['Theta2x'], acceleration_df['Alpha4'], label='Alpha4')
plt.plot(acceleration_df['Theta2x'], acceleration_df['Alpha5'], label='Alpha5')
plt.xlabel('Theta2x (degrees)')
plt.ylabel('Angular Accelerations (rad/s^2)')
plt.title('Acceleration Analysis')
plt.legend()

# Plot for Acceleration of Center of Masses
plt.subplot(3, 1, 3)
plt.plot(acceleration_cm_df['Theta2x'], acceleration_cm_df['aG2x'], label='aG2x')
plt.plot(acceleration_cm_df['Theta2x'], acceleration_cm_df['aG2y'], label='aG2y')
plt.plot(acceleration_cm_df['Theta2x'], acceleration_cm_df['aG3x'], label='aG3x')
plt.plot(acceleration_cm_df['Theta2x'], acceleration_cm_df['aG3y'], label='aG3y')
plt.plot(acceleration_cm_df['Theta2x'], acceleration_cm_df['aG4x'], label='aG4x')
plt.plot(acceleration_cm_df['Theta2x'], acceleration_cm_df['aG4y'], label='aG4y')
plt.plot(acceleration_cm_df['Theta2x'], acceleration_cm_df['aG5x'], label='aG5x')
plt.plot(acceleration_cm_df['Theta2x'], acceleration_cm_df['aG5y'], label='aG5y')
plt.plot(acceleration_cm_df['Theta2x'], acceleration_cm_df['aG6x'], label='aG6x')
plt.plot(acceleration_cm_df['Theta2x'], acceleration_cm_df['aG6y'], label='aG6y')
plt.xlabel('Theta2x (degrees)')
plt.ylabel('Acceleration (m/s^2)')
plt.title('Acceleration of Center of Masses')
plt.legend()

# Setup the figure and subplots
fig, axs = plt.subplots(4, 1, figsize=(8, 12))  # 4 subplots for each magnitude, all in one column

# Plot each magnitude in a separate subplot
axs[0].plot(theta2x_list, FO2_list, label='F02', color='blue')
axs[0].set_title('Magnitude of F02 vs Theta2x')
axs[0].set_ylabel('F02 Magnitude')

axs[1].plot(theta2x_list, FA_list, label='FA', color='green')
axs[1].set_title('Magnitude of FA vs Theta2x')
axs[1].set_ylabel('FA Magnitude')

axs[2].plot(theta2x_list, FB_list, label='FB', color='red')
axs[2].set_title('Magnitude of FB vs Theta2x')
axs[2].set_ylabel('FB Magnitude')

axs[3].plot(theta2x_list, FO4_list, label='F04', color='cyan')
axs[3].set_title('Magnitude of F04 vs Theta2x')
axs[3].set_ylabel('F04 Magnitude')

# Add legends and adjust layout
for ax in axs:
    ax.legend()
    ax.grid(True)

# Setup the figure and subplots
fig, axs = plt.subplots(3, 1, figsize=(8, 12))  # 3 subplots for each magnitude, all in one column

axs[0].plot(theta2x_list, FC_list, label='FC', color='magenta')
axs[0].set_title('Magnitude of FC vs Theta2x')
axs[0].set_ylabel('FC Magnitude')

axs[1].plot(theta2x_list, FD_list, label='FD', color='yellow')
axs[1].set_title('Magnitude of FD vs Theta2x')
axs[1].set_ylabel('FD Magnitude')

axs[2].plot(theta2x_list, T12_list, label='T12', color='orange')
axs[2].set_title('Magnitude of T12 vs Theta2x')
axs[2].set_ylabel('T12 Magnitude')

# Add legends and adjust layout
for ax in axs:
    ax.legend()
    ax.grid(True)

# Automatically adjust the layout of the figure
plt.tight_layout()

# Display the figure
plt.show()
