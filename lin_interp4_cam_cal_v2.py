# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 08:17:21 2022

Interpolation of the Navigation states at two given times

@author: mhansen
"""
import numpy as np
import math

class nothing: pass

def Q_INIT():
    q = nothing()
    q.s = 0.0
    q.v = np.array((0.0, 0.0, 0.0))

    return q

def QxQ (q1, q2):

    q_out = Q_INIT()

    q_out.s = q1.s*q2.s - (q1.v[0]*q2.v[0] + q1.v[1]*q2.v[1] + q1.v[2]*q2.v[2])
    q_out.v[0] = q1.s*q2.v[0] + q2.s*q1.v[0] - (q1.v[1]*q2.v[2] - q1.v[2]*q2.v[1])
    q_out.v[1] = q1.s*q2.v[1] + q2.s*q1.v[1] - (q1.v[2]*q2.v[0] - q1.v[0]*q2.v[2])
    q_out.v[2] = q1.s*q2.v[2] + q2.s*q1.v[2] - (q1.v[0]*q2.v[1] - q1.v[1]*q2.v[0])

    return q_out

def Q_CONJ(q_in):

    q_out = Q_INIT()

    q_out.s = q_in.s
    q_out.v[0] = -q_in.v[0]
    q_out.v[1] = -q_in.v[1]
    q_out.v[2] = -q_in.v[2]
    
    return q_out

def Q_NORM(q):
    q_mag = math.sqrt(q.s*q.s + q.v[0]*q.v[0] + q.v[1]*q.v[1] + q.v[2]*q.v[2])
    if( q_mag > 0.01 ):
        q_mag = 1.0/q_mag
    else:
        q_mag = 0.0
        
    q.s *= q_mag 
    q.v[0] *= q_mag 
    q.v[1] *= q_mag 
    q.v[2] *= q_mag 

    return q

def Q_to_DCM(Q):
    """
    Transformation matrix from quaternion (T = I - 2*s*[v x] + [v x]^2)
    """
    s = Q.s; v0 = Q.v[0]; v1 = Q.v[1]; v2= Q.v[2]
    
    T = np.mat(np.zeros((3,3)))
    T[0,0] = 1.0 - 2.0*(v2*v2 + v1*v1)
    T[0,1] = 2.0*(v1*v0 + s*v2)
    T[0,2] = 2.0*(v2*v0 - s*v1)#
    T[1,0] = 2.0*(v1*v0 - s*v2)
    T[1,1] = 1.0 - 2.0*(v2*v2 + v0*v0)
    T[1,2] = 2.0*(v2*v1 + s*v0)
    T[2,0] = 2.0*(v2*v0 + s*v1)
    T[2,1] = 2.0*(v2*v1 - s*v0)
    T[2,2] = 1.0 - 2.0*(v1*v1 + v0*v0)
        
    return T

def mat2eul123(dcm):
    '''%
    %/* FUNCTION:    A transformation matrix is generated using the input
    %        ROLL-PITCH-YAW euler angles (in radians). The angles
    %        are input as:
    %        angle[0] = ROLL     - M_PI   <= ROLL  < M_PI
    %        angle(2) = PITCH    - M_PI/2 <= PITCH < M_PI/2
    %        angle(3) = YAW      - M_PI   <= YAW   < M_PI
    %*/
    '''

    s1 = 0;		#/* SINE OF ROLL */
    c1 = 0; 		#/* COSINE OF ROLL */
    s2 = 0;		#/* SINE OF PITCH */
    c2 = 0;		#/* COSINE OF PITCH */
    s3 = 0;		#/* SINE OF YAW */
    c3 = 0;		#/* COSINE OF YAW */

    angle = np.array((0.0,0.0,0.0))

    angle[1] = np.arcsin( dcm[2,0] ) ;
    
    
    if( abs( angle[1] - np.pi/2.0 ) < 1.0e-6 ):

        print('WARNING: Singularity in deuler_123() at PITCH = 90.0 deg')

        angle[0] = np.arctan2( dcm[0,1], dcm[1,1] ) ;
        angle[2] = 0.0 ;

    elif( abs( angle[1] + np.pi/2.0 ) < 1.0e-6 ):
        print('WARNING: Singularity in deuler_123() at PITCH = -90.0 deg')

        angle[0] = np.arctan2( -dcm[0,1], dcm[1,1] ) ;
        angle[2] = 0.0 ;
  
    else:

        angle[0] = np.arctan2( -dcm[2,1], dcm[2,2] ) ;
        angle[2] = np.arctan2( -dcm[1,0], dcm[0,0] ) ;


    return angle

# Input Parameters
case2inrtl  = 0  #if 0 assume quaternion is inrtl to case
chassis2opt = 0  #if 1 transform quaternion to an inrtl to opt frame
theta_cant  = 10*np.pi/180 #(deg)

q_case1_to_inrtl = Q_INIT()
q_case2_to_inrtl = Q_INIT()
q_rotate = Q_INIT()
q_temp   = Q_INIT()

q_rotate.s    = np.cos(theta_cant/2)
q_rotate.v[0] = np.sin(theta_cant/2)
q_rotate.v[1] = 0.0
q_rotate.v[2] = 0.0

q_case1_to_inrtl.s =   0.786066  #-0.785230 
q_case1_to_inrtl.v[0] = 0.531661 # 0.532809
q_case1_to_inrtl.v[1] = 0.278631 # 0.276865
q_case1_to_inrtl.v[2] =-0.147651 #-0.151242
print(" ")
print("----------- Quaternion Rotation Analysis -----------")
print("Expected Rotation: ", 90, " deg")
print("qr_norm_1: ", np.sqrt(q_case1_to_inrtl.s*q_case1_to_inrtl.s + 
                           q_case1_to_inrtl.v[0]*q_case1_to_inrtl.v[0] + 
                           q_case1_to_inrtl.v[1]*q_case1_to_inrtl.v[1] +
                           q_case1_to_inrtl.v[2]*q_case1_to_inrtl.v[2] ))
if chassis2opt:
    q_temp.s = q_case1_to_inrtl.s
    q_temp.v[0] = q_case1_to_inrtl.v[0]
    q_temp.v[1] = q_case1_to_inrtl.v[1]
    q_temp.v[2] = q_case1_to_inrtl.v[2]
    
    if case2inrtl:
        q_rotate = Q_CONJ(q_rotate)
        q_case1_to_inrtl = QxQ(q_temp, q_rotate)
    else:
        q_case1_to_inrtl = QxQ(q_rotate, q_temp)
        
Q_NORM(q_case1_to_inrtl)

q_case2_to_inrtl.s =    0.535216  #-0.400954
q_case2_to_inrtl.v[0] = 0.602515  # 0.139360
q_case2_to_inrtl.v[1] = 0.006196  # 0.590585
q_case2_to_inrtl.v[2] = 0.592014  #-0.686311

if chassis2opt:
    q_temp.s = q_case2_to_inrtl.s
    q_temp.v[0] = q_case2_to_inrtl.v[0]
    q_temp.v[1] = q_case2_to_inrtl.v[1]
    q_temp.v[2] = q_case2_to_inrtl.v[2]
    
    if case2inrtl:
        q_rotate = Q_CONJ(q_rotate)
        q_case2_to_inrtl = QxQ(q_temp, q_rotate)
    else:
        q_case2_to_inrtl = QxQ(q_rotate, q_temp)
    
print("qr_norm_2: ", np.sqrt(q_case2_to_inrtl.s*q_case2_to_inrtl.s + 
                           q_case2_to_inrtl.v[0]*q_case2_to_inrtl.v[0] + 
                           q_case2_to_inrtl.v[1]*q_case2_to_inrtl.v[1] +
                           q_case2_to_inrtl.v[2]*q_case2_to_inrtl.v[2] ))
Q_NORM(q_case2_to_inrtl)

# Quaternion Multiply
if case2inrtl:
    q_inrtl_to_case2 = Q_CONJ(q_case2_to_inrtl)
    q_case1_to_case2 = QxQ(q_inrtl_to_case2, q_case1_to_inrtl)
else:
    q_inrtl_to_case1 = Q_CONJ(q_case1_to_inrtl)
    q_case1_to_case2 = QxQ(q_case2_to_inrtl, q_inrtl_to_case1)
    
T_case1_to_case2 = Q_to_DCM(q_case1_to_case2)
eul_case1_2_case2 = mat2eul123(T_case1_to_case2)

# Determine Rotaion magnitude as determined by ST
mag_quat = 2*np.arccos(q_case1_to_case2.s)*180/np.pi
q_vec_norm = np.sqrt(q_case1_to_case2.v[0]*q_case1_to_case2.v[0] + 
                  q_case1_to_case2.v[1]*q_case1_to_case2.v[1]+
                  q_case1_to_case2.v[2]*q_case1_to_case2.v[2])

print("Rotation Magnitude by Quaternion: ", mag_quat)
print(" ")
print("Quat case1 to case2: ", q_case1_to_case2.s, q_case1_to_case2.v)
print("Eigen Vector Between Quaternion: ", q_case1_to_case2.v/q_vec_norm)
print("Euler Angle RPY: ", eul_case1_2_case2*180/np.pi)