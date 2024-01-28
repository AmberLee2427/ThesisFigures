"""
Creating the rotation curve of the Milky Way based on distances, 
proper motions, and radial velocities of Cepheids.

Usage: 

P. Mroz @ OAUW, 12 Sep 2018

"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import time

def calculate_galactocentric_velocity_and_radius (ra,dec,l,b,dist,\
    dist_err,pm_ra,pm_ra_err,pm_dec,pm_dec_err,pm_corr,v_rad,v_rad_err,\
    ra_G,dec_G,R_0,U_0,V_0,W_0,THETA_0):
    """
    Calculating Galactocentric velocity and radius
    
    Parameters
    ----------
    ra : right ascension (radians)
    dec : declination (radians)
    l : Galactic longitude (radians)
    b : Galactic latitude (radians)
    dist : Heliocentric distance (kpc)
    dist_err : uncertainty of Heliocentric distance (kpc)
    pm_ra : proper motion in right ascension (mas/yr)
    pm_ra_err : uncertainty of proper motion in right ascension (mas/yr)
    pm_dec : proper motion in declination (mas/yr)
    pm_dec_err : uncertainty of proper motion in declination (mas/yr)
    pm_corr : correlation coefficient between pm_ra and pm_dec
    v_rad : Heliocentric radial velocity (km/s)
    v_rad_err : uncertainty of Heliocentric radial velocity (km/s)
    ra_G : right ascension of the North Galactic pole (radians)
    dec_G : declination of the North Galactic pole (radians)
    R_0 : distance to the Galactic center (kpc)
    U_0 : solar motion toward the Galactic center (km/s)
    V_0 : solar motion toward l=90 deg (km/s)
    W_0 : solar motion toward the North Galactic pole (km/s)
    THETA_0 : rotation speed of the LSR (km/s)
    
    Returns
    -------
    R : Galactocentric distance (kpc)
    R_err : uncertainty of Galactocentric distance (kpc) 
    V_TOT : circular velocity relative to the GC (km/s) 
    V_TOT_ERR : uncertainty of the circular velocity (km/s) 
    """
    # Transforming equatorial proper motions to the Galactic system
    C1 = np.sin(dec_G)*np.cos(dec)-np.cos(dec_G)*np.sin(dec)*np.cos(ra-ra_G)
    C2 = np.cos(dec_G)*np.sin(ra-ra_G)
    norm = np.sqrt(C1*C1+C2*C2)
    C1 /= norm
    C2 /= norm
    pm_l = C1*pm_ra+C2*pm_dec
    pm_b = C1*pm_dec-C2*pm_ra
    # Transforming error bars of proper motions
    J1 = np.array([[C1,C2],[-C2,C1]])
    J2 = np.transpose(J1)
    CORR = np.array([[pm_ra_err**2,pm_corr*pm_ra_err*pm_dec_err], \
    [pm_corr*pm_ra_err*pm_dec_err,pm_dec_err**2]])
    COV = np.dot(np.dot(J1,CORR),J2)
    pm_l_err = np.sqrt(COV[0,0])
    pm_b_err = np.sqrt(COV[1,1])
    pm_lb_corr = COV[0,1] / pm_l_err / pm_b_err
    # Calculating tangential components of velocity vector and their error bars
    v_l = 4.74 * dist * pm_l
    v_b = 4.74 * dist * pm_b
    v_l_err = np.absolute(v_l) * np.sqrt((dist_err/dist)**2 + (pm_l_err/pm_l)**2)
    v_b_err = np.absolute(v_b) * np.sqrt((dist_err/dist)**2 + (pm_b_err/pm_b)**2)
    # Calculating the velocity vector in the Cartesian Galactic coordinates
    # at the location of the Sun:
    C11 = np.cos(b)*np.cos(l)
    C12 = (-1.0)*np.sin(b)*np.cos(l)
    C13 = (-1.0)*np.sin(l)
    C21 = np.cos(b)*np.sin(l)
    C22 = (-1.0)*np.sin(b)*np.sin(l)
    C23 = np.cos(l)
    C31 = np.sin(b)
    C32 = np.cos(b)
    C33 = 0.0
    U_1 = C11*v_rad + C12*v_b + C13*v_l
    V_1 = C21*v_rad + C22*v_b + C23*v_l
    W_1 = C31*v_rad + C32*v_b + C33*v_l
    # Adding orbital motion of the Sun 
    U_2 = U_1 + U_0
    V_2 = V_1 + V_0 + THETA_0
    W_2 = W_1 + W_0
    # Calculating the covariance matrix
    J1 = np.array([[C11,C12,C13],[C21,C22,C23],[C31,C32,C33]])
    J2 = np.transpose(J1)
    CORR = np.array([[v_rad_err**2,0.0,0.0],[0.0,v_b_err**2,\
    pm_lb_corr*v_b*v_l*(pm_l_err/pm_l)*(pm_b_err/pm_b)],\
    [0.0,pm_lb_corr*v_b*v_l*(pm_l_err/pm_l)*(pm_b_err/pm_b),v_l_err**2]])
    TMP = np.dot(CORR,J2)
    COV = np.dot(J1,TMP)
    U_2_err = np.sqrt(COV[0,0])
    V_2_err = np.sqrt(COV[1,1])
    W_2_err = np.sqrt(COV[2,2])
    # Calculating total velocity and its error bar
    V_TOT = np.sqrt(U_2**2 + V_2**2 + W_2**2)
    V_TOT_ERR = U_2*U_2*COV[0,0] + V_2*V_2*COV[1,1] + W_2*W_2*COV[2,2] \
    + 2.0*U_2*V_2*COV[1,0] + 2.0*U_2*W_2*COV[2,0] + 2.0*V_2*W_2*COV[1,2]
    V_TOT_ERR = np.sqrt(V_TOT_ERR) / V_TOT
    # Calculating Galactocentric distance
    R = np.sqrt(R_0*R_0 + (dist*np.cos(b))**2 - 2.0*R_0*dist*np.cos(b)*np.cos(l))
    R_err = (dist_err / R)*np.cos(b)*np.absolute(dist*np.cos(b)-R_0*np.cos(l))
    # Calculate angle beta between the Sun and the source as viewed from the GC
    sinbeta = dist*np.cos(b)*np.sin(l)/R
    cosbeta = (R_0-dist*np.cos(b)*np.cos(l))/R
    beta = np.arctan2(sinbeta,cosbeta)
    U_s = U_2*np.cos(beta)-V_2*np.sin(beta)
    V_s = V_2*np.cos(beta)+U_2*np.sin(beta)
    W_s = W_2
    J1 = np.array([[np.cos(beta),-np.sin(beta),0.0],[np.sin(beta),np.cos(beta),0.0],[0.0,0.0,1.0]])
    J2 = np.transpose(J1)
    CORR = np.dot(J1,np.dot(COV,J2))
    U_s_err = np.sqrt(CORR[0,0])
    V_s_err = np.sqrt(CORR[1,1])
    W_s_err = np.sqrt(CORR[2,2])
    return R, R_err, V_s, V_s_err
    

def main():

    matplotlib.rcParams['text.latex.preamble'] = [
           r'\usepackage{helvet}',    # set the normal font here
           r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
           r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
    ]  

    plt.rc('pdf', fonttype=42)
    plt.rc('text', usetex=True)


    #-----Parameters and constants-----------------------------------------#

    deg = np.pi / 180.0 # degrees to radians
    dec_G = 27.12825 * deg # North Galactic pole
    ra_G = 192.85948 * deg # North Galactic pole

    #R0 = 8.122 # distance to the Galactic center, Abuter et al. 2018 (Gravity)
    # velocity of the Sun relative to the LSR, Schonrich, Binney & Dehnen 2010:
    #U0 = 11.1
    #V0 = 12.24
    #W0 = 7.25
    #VEL_0 = 230.0 # rotation speed of LSR

    # MODEL 2:
    R0 = 8.09 # distance to the Galactic center
    # velocity of the Sun relative to the LSR:
    U0 = 10.1 
    V0 = 12.3
    W0 = 7.3
    VEL_0 = 233.6

    #-----Reading data from file-------------------------------------------#

    filename = 'Cepheids.dat'
    dtype = 'S20,int,float,float,float,float,float,float,float,float, float,float,float,float,float'
    try:
        ident, flag, ra_ogle, dec_ogle, glon_ogle, glat_ogle, dist, \
        dist_err, pm_ra, pm_ra_err, pm_dec, pm_dec_err, pm_corr, v_r, \
        v_r_err = np.loadtxt(filename, unpack=True, dtype=dtype)
        N = len(ra_ogle)
    except IOError:
        exit('File %s not found. Exiting.'%filename)
        
    l = glon_ogle * deg # Galactic longitude
    b = glat_ogle * deg # Galactic latitude
    ra = ra_ogle * deg # right ascension
    dec = dec_ogle * deg # declination

    #-----Calculating Galactocentric velocity------------------------------#

    R, R_err, V_TOT, V_TOT_ERR = [], [], [], []

    for i in xrange(N):
        _R,_R_err,_V_TOT,_V_TOT_ERR = calculate_galactocentric_velocity_and_radius \
        (ra[i],dec[i],l[i],b[i],dist[i],dist_err[i],pm_ra[i],pm_ra_err[i], \
        pm_dec[i],pm_dec_err[i],pm_corr[i],v_r[i],v_r_err[i],ra_G,dec_G,R0,U0,V0,W0,VEL_0)
        R.append(_R)
        R_err.append(_R_err)
        V_TOT.append(_V_TOT)
        V_TOT_ERR.append(_V_TOT_ERR)

    R = np.array(R)
    R_err = np.array(R_err)
    V_TOT = np.array(V_TOT)
    V_TOT_ERR = np.array(V_TOT_ERR)

    #-----Saving data to file----------------------------------------------#

    fname = 'rotation_curve.txt'
    with open(fname, 'w') as f:
        f.write('# Galactic rotation curve based on Classical Cepheids\n')
        f.write('# (Mroz et al. 2018)\n')
        f.write('# Using R_0 = %.2f kpc, Omega_0 = %.1f km/s,\n'%(R0,VEL_0))
        f.write('# U_0 = %.1f km/s, V_0 = %.1f km/s, W_0 = %.1f km/s\n'%(U0,V0,W0))
        str_time = time.strftime('%b %d %Y %H:%M %Z')
        f.write('# Created: %s\n'%str_time)
        f.write('# Name R R_err V V_err\n')
        for i in xrange(N):
            if flag[i] == 0: continue
            f.write('%-20s %6.3f %5.3f %6.2f %5.2f\n'%(ident[i],R[i],R_err[i],V_TOT[i],V_TOT_ERR[i]))

    #------Creating rotation curve----------------------------------------#

    fig = plt.figure(figsize=(7,4))
    ax = plt.gca()
    filtr = (flag == 1)
    ax.errorbar(R[filtr], V_TOT[filtr], yerr=V_TOT_ERR[filtr], xerr=R_err[filtr], fmt='o', ms=1.5, elinewidth=0.5, color='navy', zorder=2)
        
    ax.set_xlim(0.0, 20.0)
    ax.set_ylim(0.0, 350.0)
    ax.tick_params(axis='both',which='both',direction='in',top=True,right=True)
    ax.set_xticks(np.arange(0.0,20.0,1.0),minor=True)
    ax.set_xticks(np.arange(0.0,21.0,2.0),minor=False)
    ax.set_yticks(np.arange(0.0,350.0,10.0),minor=True)
    ax.set_xlabel(r'Distance from the Galactic center (kpc)', fontsize=12)
    ax.set_ylabel(r'Circular velocity (km/s)', fontsize=12)

    #ax.scatter((13.5),(30.0),s=4,color='navy')
    #ax.text(14.0, 30.0, r'This work', fontsize=12, color='navy',verticalalignment='center')

    plt.savefig('rotation_curve.pdf',bbox='tight',dpi=200)
    
if __name__ == '__main__':
    main()