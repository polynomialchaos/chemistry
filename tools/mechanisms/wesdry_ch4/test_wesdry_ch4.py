import cantera as ct

gas = ct.Solution( 'wesdry_ch4.cti' )
gas.TPY = 1600.0, 1e5, {'CH4':0.1, 'O2':0.1, 'Ar':0.8}

gas()