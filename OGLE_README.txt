=======================================================================

                  Wide-Orbit Exoplanets are Common.
     Analysis of Nearly 20 Years of OGLE Microlensing Survey Data

            R. Poleski, J. Skowron, P. Mróz, A. Udalski,
    M.K. Szymanski, P. Pietrukowicz, K. Ulaczyk, K. Rybicki,
                 P. Iwanek, M.Wrona and M. Gromadzki

          Acta Astronomica 71, 1 (2021); arXiv: 2104.02079

=======================================================================
rpoleski@astrouw.edu.pl

This directory presents results of analysis of occurrence of wide-orbit
microlensing planets.

The directory structure is as follows:

README            - this file,

fields.dat        - information on fields studied,

events.dat        - list of microlensing events analyzed,

detection_eff.dat - detection efficiency

post_def.dat      - posterior samples for default selection criteria,

post_ext.dat      - posterior samples for extended selection criteria,

post_def_add_ob170114.dat - posterior samples for default selection
                            criteria and assuming wide-orbit solution for
                            OGLE-2017-BLG-0114,

post_ext_add_ob170114.dat - posterior samples for extended selection
                            criteria and assuming wide-orbit solution for
                            OGLE-2017-BLG-0114.


================================================================================
Byte-by-byte Description of file: fields.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label      Explanations
--------------------------------------------------------------------------------
   1-  6 A6     ---     Name       Name of the field
   8-  9 I2     h       RAh        Right Ascension J2000 (hours)
  11- 12 I2     min     RAm        Right Ascension J2000 (minutes)
  14- 17 F4.1   s       RAs        Right Ascension J2000 (seconds)
  19- 21 I3     deg     DEd        Declination J2000 (degrees)
  23- 24 I2     arcmin  DEm        Declination J2000 (arcminutes)
  26- 27 I2     arcsec  DEs        Declination J2000 (arcseconds)
      29 A1     ---     OGLEphase  OGLE phase: 3 or 4
  31- 37 F7.2   d       DeltaTime  Total time of observations
  39- 43 I5     ---     Nepochs    Number of epochs
  45- 50 F6.4   d       MedianDT   Median time difference between consecutive epochs
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

================================================================================
Byte-by-byte Description of file: events.dat
--------------------------------------------------------------------------------
   Bytes Format Units    Label      Explanations
--------------------------------------------------------------------------------
   1- 16 A16    ---      Name       Event name
  18- 19 I2     h        RAh        Right Ascension J2000 (hours)
  21- 22 I2     min      RAm        Right Ascension J2000 (minutes)
  24- 28 F5.2   s        RAs        Right Ascension J2000 (seconds)
  30- 32 I3     deg      DEd        Declination J2000 (degrees)
  34- 35 I2     arcmin   DEm        Declination J2000 (arcminutes)
  37- 40 F4.1   arcsec   DEs        Declination J2000 (arcseconds)
  42- 42 I1     ---      OGLEphase  Phase of OGLE project: 3 or 4
  44- 61 A18    ---      EWSid      ID on EWS website or "-" if not announced there
  63- 67 I5     ---      Nepochs    Number of epochs
  69- 76 F8.5   ---      chi2dof    Chi^2 per number of degrees of freedom
  78- 90 F13.5  ---      t0         HJD of minimum separation
  92-100 F9.6   ---      u0         Impact parameter relative to theta_E
 102-111 F10.6  d        tE         Einstein ring crossing time
 113-120 F8.5   ---      piEN       North component of microlensing parallax
 122-129 F8.5   ---      piEE       East component of microlensing parallax
 131-137 I7     ---      t0par      HJD of parameter reference epoch
 139-144 F6.3   mag      Ibase      Baseline brightness in I band
 146-151 F6.3   mag      Isource    Source brightness in I band
 153-158 F6.3   mag      IRC        Red clump brightness in I band
 160-160 A1     ---    f_IRC        Flag for IRC (1)
 162-169 F8.6   ---      rho        Source radius relative to theta_E
 171-178 F8.6   ---    e_rho        Error on rho (lower value)
 180-187 F8.6   ---    E_rho        Error on rho (upper value)
 189-189 A1     ---    f_rho        Flag for rho (2)
--------------------------------------------------------------------------------
Note (1):
    N = value from Nataf et al. (2013; ApJ 769, 88)
    f = value fitted by the authors
    - = no data
Note (2):
    m = value from Galactic model
    f = value from microlensing fit
    - = no data
--------------------------------------------------------------------------------

================================================================================
Byte-by-byte Description of file: detection_eff.dat
--------------------------------------------------------------------------------
   Bytes Format Units Label               Explanations
--------------------------------------------------------------------------------
   1-  3 F3.1   ---     s                   Projected separation relative to Einstein ring radius
   5- 12 F8.6   ---     q                   Mass-ratio
  14- 20 F7.3   ---     SensitivityDefault  Survey sensitivity for default detection criteria (Fig. 6 of the paper)
  22- 28 F7.3   ---     SensitivityExtend   Survey sensitivity for extended detection criteria
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

================================================================================
Byte-by-byte Description of files: post_def.dat, post_ext.dat, 
post_def_add_ob170114.dat, and post_ext_add_ob170114.dat
--------------------------------------------------------------------------------
   1-  8 F8.6   ---     A             Multiplicative constant
  10- 18 F9.6   ---     n             Exponent of s
  20- 28 F9.6   ---     m             Exponent of q
  30- 38 F9.6   ---     log10(prior)  Log10 of prior probability
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

Any presentation of the scientific analysis or usage of the data from this
directory should include appropriate reference to the paper.
