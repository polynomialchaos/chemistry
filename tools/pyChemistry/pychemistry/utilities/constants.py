################################################################################
# @file constants.py
# @author Florian Eigentler
# @brief
# @version 1.0.0
# @date 2022-02-22
# @copyright Copyright (c) 2022 by Florian Eigentler.
#  This work is licensed under terms of the MIT license (<LICENSE>).
################################################################################
NA = 6.022140857e23  # Avogadro number (mol-1)
ENH_VW = 2.2         # Van Waals enhancement factor (-)

KB = 1.38064852e-23  # Boltzmann constant (m2 kg s-2 K-1)
RM = NA * KB         # Universal gas constant (mol-1 m2 kg s-2 K-1)
P0 = 101325.0        # Standard pressure (Pa)
T0 = 273.15          # Standard temperature (K)

CM_M = 0.01                 # Conversion from (cm) to (m)
CAL_JOULE = 4.184           # Conversion from (cal) to (J)
JOULE_KELVINS = RM          # Conversion from (J) to (K)
ANGSTROM_SI = 10.0**(-10)   # Conversion from (Angstrom) to (m)
DEBYE_SI = 10.0**(-24.5)    # Conversion from (Debye) to (m3 m-2 J J-2)

REF_ELEMENTS = {
    'H':  0.00100794,
    'HE':  0.0040026022,
    'C':  0.0120107,
    'N':  0.0140067,
    'O':  0.0159994,
    'NE':  0.0201797,
    'AR':  0.039948,
}
