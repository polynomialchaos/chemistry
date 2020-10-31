#!/usr/bin/env python3
####################################################################################################################################
# Chemistry - Finite Rate Chemistry library
# (c) 2020 | Florian Eigentler
####################################################################################################################################
import argparse, logging, os, sys

from pychemistry import parse_input, write_chemistry_data

####################################################################################################################################
# Definitions
# ----------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Class Definitions
# ----------------------------------------------------------------------------------------------------------------------------------
class LessThanFilter(logging.Filter):
    def __init__( self, exclusive_maximum, name='' ):
        super( LessThanFilter, self ).__init__( name )
        self.max_level = exclusive_maximum

    def filter( self, record ):
        #non-zero return means we log this message
        return 1 if record.levelno < self.max_level else 0

logging.getLogger().setLevel( logging.DEBUG )

err_handler = logging.StreamHandler( stream=sys.stderr )
err_handler.setLevel( logging.ERROR )
err_handler.setFormatter( logging.Formatter( '%(levelname)s:%(name)s:%(message)s' ) )
logging.getLogger().addHandler( err_handler )

out_handler = logging.StreamHandler( stream=sys.stdout )
out_handler.setLevel( logging.INFO )
out_handler.addFilter( LessThanFilter( logging.ERROR ) )
out_handler.setFormatter( logging.Formatter( '%(levelname)s:%(name)s:%(message)s' ) )
logging.getLogger().addHandler( out_handler )

####################################################################################################################################
# Functions
# ----------------------------------------------------------------------------------------------------------------------------------
def main():

    # define the argument parser
    parser = argparse.ArgumentParser( description='Chemistry - Finite Rate Chemistry library - preprocessing' )
    parser.add_argument( '-d', '--debugging', dest='debugging', action='store_true', help='Debugging output' )
    parser.add_argument( '-i', '--input', dest='input', type=str, required=True, help='The chemkin fortmat input file' )
    parser.add_argument( '-th', '--thermo', dest='thermo', type=str, default=None, help='The chemkin fortmat thermo file' )
    parser.add_argument( '-tr', '--transport', dest='transport', default=None, type=str, help='The chemkin fortmat transport file' )
    parser.add_argument( '-o', '--output', dest='output', type=str, help='The Chemistry output file' )
    parser.add_argument( '-ck', '--chemkin', dest='chemkin', type=str, help='The chemkin format output file' )
    args = parser.parse_args()

    # setup additional output for debugging
    if args.debugging:
        debug_handler = logging.FileHandler( '{:}.log'.format( args.inifile ), mode='w' )
        debug_handler.setLevel( logging.DEBUG )
        debug_handler.setFormatter( logging.Formatter( '%(levelname)s:%(name)s:%(message)s' ) )
        logging.getLogger().addHandler( debug_handler )

    # read the mechanism
    mechanism = parse_input( mech_path=args.input, thermo_path=args.thermo, transport_path=args.transport )

    # check the parsed data from the input file
    if mechanism.is_valid(): logging.info( 'Mechanism passed verification!' )
    else: logging.warning( 'Mechanism failed verification!' )

    # write the Chemistry format
    output = args.output if args.output else '{:}.h5'.format( os.path.basename( args.input ) )
    write_chemistry_data( output, mechanism )

    # write the chemkin format
    if args.chemkin:
        mechanism.chemkinify( args.chemkin, unit_Ea='KELVINS' )

####################################################################################################################################
# CALL BY SCRIPT
# ----------------------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    main()