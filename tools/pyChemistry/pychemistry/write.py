####################################################################################################################################
# pyChemistry - Python package for FV3D preprocessing
# (c) 2020 | Florian Eigentler
####################################################################################################################################
import logging, h5py, datetime
import numpy as np

####################################################################################################################################
# Definitions
# ----------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Class Definitions
# ----------------------------------------------------------------------------------------------------------------------------------

####################################################################################################################################
# Functions
# ----------------------------------------------------------------------------------------------------------------------------------
def get_specii_reaction_points( specii, reactions ):
    """Generate the list of reaction points."""

    reac_points, reac_nus = [], []
    for species in specii:
        sp_reac_points, sp_reac_nus = [], []
        for idx_re, reaction in enumerate( reactions ):
            found_species, reac_nu = False, 0.0

            for key in reaction.reactants:
                if key == species:
                    found_species   = True
                    reac_nu        -= reaction.reactants[key]

            for key in reaction.products:
                if key == species:
                    found_species   = True
                    reac_nu        += reaction.products[key]

            if found_species:
                sp_reac_points.append( idx_re )
                sp_reac_nus.append( reac_nu )

        reac_points.append( sp_reac_points )
        reac_nus.append( sp_reac_nus )

    return reac_points, reac_nus

def write_chemistry_data( path, mechanism ):
    """Generate the Chemistry format mechanism file."""
    logging.info( 'Write Chemistry format mechanism file' )

    # --- elements ---
    elSymbol                = [x for x in mechanism.elements]
    elAtomicMass            = [mechanism.elements[x].mass for x in elSymbol]

    # --- species ---
    spSymbol                = [x for x in mechanism.specii]

    tmpInerts               = mechanism.inert_specii()
    spIsInert               = [x in tmpInerts for x in spSymbol]

    sortIdx                 = [x[0] for x in sorted( enumerate( spIsInert ), key=(lambda x: x[1] is True) )]
    spSymbol                = [spSymbol[x] for x in sortIdx]
    spIsInert               = [spIsInert[x] for x in sortIdx]
    spSymbolM               = spSymbol + ['M']

    # --- thermo ---
    thDatas                 = [mechanism.thermos[symbol] for symbol in spSymbol]
    getElCount              = (lambda x: [x[symbol] if symbol in x else 0.0 for symbol in elSymbol])
    thComposition           = [getElCount( dict( x.composition ) ) for x in thDatas]
    thPhase                 = [x.phase for x in thDatas]
    thBounds                = [x.bounds for x in thDatas]
    thCoeffHigh             = [x.coeff_high for x in thDatas]
    thCoeffLow              = [x.coeff_low for x in thDatas]

    # --- transport ---
    trDatas                 = [mechanism.transports[symbol] for symbol in spSymbol]
    trGeom                  = [x.geom for x in trDatas]
    trPotLJ                 = [x.pot_lj for x in trDatas]
    trColLJ                 = [x.col_lj for x in trDatas]
    trDipMo                 = [x.dip_mo for x in trDatas]
    trPol                   = [x.pol for x in trDatas]
    trRotRel                = [x.rot_rel for x in trDatas]

    # --- species-reaction points ---
    spReacPoints, spReacNus = get_specii_reaction_points( spSymbol, mechanism.reactions )
    spNReacPoints           = [len( x ) for x in spReacPoints]

    # --- reactions ---
    reType                  = [x.reaction_type.value for x in mechanism.reactions]
    reFallOffSpecies        = [spSymbolM.index( x.falloff_species ) if x.falloff_species else -1 for x in mechanism.reactions]

    reArrCoeff              = [x.arr_coeff for x in mechanism.reactions]
    reIsReversible          = [x.is_reversible for x in mechanism.reactions]

    reNReactants            = [len( x.reactants ) for x in mechanism.reactions]
    reReactants             = [[spSymbolM.index( y ) for y in x.reactants] for x in mechanism.reactions]
    reNuReactants           = [[x.reactants[y] for y in x.reactants] for x in mechanism.reactions]
    reOrdReactants          = [x.reactant_orders for x in mechanism.reactions]

    reNProducts             = [len( x.products ) for x in mechanism.reactions]
    reProducts              = [[spSymbolM.index( y ) for y in x.products] for x in mechanism.reactions]
    reNuProducts            = [[x.products[y] for y in x.products] for x in mechanism.reactions]
    reOrdProducts           = [x.product_orders for x in mechanism.reactions]

    reSumNu                 = [x.sum_nu for x in mechanism.reactions]

    reHasRevArr             = [True if x.rev_arr_coeff else False for x in mechanism.reactions]
    reRevArrCoeff           = [x.rev_arr_coeff for x in mechanism.reactions]

    reHasHighArr            = [True if x.adv_arr_key == 'HIGH' else False for x in mechanism.reactions]
    reHasLowArr             = [True if x.adv_arr_key == 'LOW' else False for x in mechanism.reactions]
    reAdvArrCoeff           = [x.adv_arr_coeff for x in mechanism.reactions]

    reHasTroe               = [True if x.troe_coeff else False for x in mechanism.reactions]
    reNTroeCoeff            = [len( x.troe_coeff ) for x in mechanism.reactions]
    reTroeCoeff             = [x.troe_coeff for x in mechanism.reactions]

    reHasEfficiencies       = [True if x.efficiencies else False  for x in mechanism.reactions]
    reNEfficiencies         = [len( x.efficiencies.keys() ) for x in mechanism.reactions]
    reSpEfficiencies        = [[spSymbolM.index( y ) for y in x.efficiencies] for x in mechanism.reactions]
    reEfficiencies          = [[x.efficiencies[y] for y in x.efficiencies] for x in mechanism.reactions]

    # --- output ---
    getMaxCount = (lambda x: max( [len( y ) for y in x] ))
    fillArray = (lambda lists, fill, maxLength=getMaxCount: [x + [fill for _ in range( maxLength( lists )-len( x ) )] for x in lists])

    with h5py.File( path, 'w' ) as h5f:
        h5f.attrs['date']   = np.string_( str( datetime.datetime.now() ) )

        # elements
        group = h5f.create_group( r'ELEMENTS' )
        group.attrs[r'n_elements'] = len( elSymbol )

        group.create_dataset( r'symbol'             , data=[np.string_( x ) for x in elSymbol] )
        group.create_dataset( r'mass'               , data=[x for x in elAtomicMass] )

        # specii
        group = h5f.create_group( r'SPECII' )
        group.attrs['n_specii'] = len( spSymbol )
        group.attrs['max_reac_points'] = max( spNReacPoints )

        group.create_dataset( r'symbol'             , data=[np.string_( x ) for x in spSymbol] )
        group.create_dataset( r'is_inert'           , data=spIsInert )

        group.create_dataset( r'composition'        , data=thComposition )
        group.create_dataset( r'phase'              , data=thPhase )
        group.create_dataset( r'bounds'             , data=thBounds )
        group.create_dataset( r'coeff_high'         , data=thCoeffHigh )
        group.create_dataset( r'coeff_low'          , data=thCoeffLow )

        group.create_dataset( r'geom'               , data=trGeom )
        group.create_dataset( r'pot_lj'             , data=trPotLJ )
        group.create_dataset( r'col_lj'             , data=trColLJ )
        group.create_dataset( r'dip_mo'             , data=trDipMo )
        group.create_dataset( r'pol'                , data=trPol )
        group.create_dataset( r'rot_rel'            , data=trRotRel )

        group.create_dataset( r'n_reac_points'      , data=spNReacPoints )
        group.create_dataset( r'reac_points'        , data=fillArray( spReacPoints, -1 ) )
        group.create_dataset( r'nu_reac_points'     , data=fillArray( spReacNus, 0.0 ) )

        # reactions
        group = h5f.create_group( r'REACTIONS' )
        group.attrs['n_reactions'] = len( reType )
        group.attrs['max_reactants'] = max( reNReactants )
        group.attrs['max_products'] = max( reNProducts )
        group.attrs['max_troe_coeff'] = max( reNTroeCoeff )
        group.attrs['max_efficiencies'] = max( reNEfficiencies )

        group.create_dataset( r'type'               , data=reType )
        group.create_dataset( r'falloff_species'    , data=reFallOffSpecies )

        group.create_dataset( r'arr_coeff'          , data=reArrCoeff )
        group.create_dataset( r'is_reversible'      , data=reIsReversible )

        group.create_dataset( r'n_reactants'        , data=reNReactants )
        group.create_dataset( r'reactants'          , data=fillArray( reReactants, -1 ) )
        group.create_dataset( r'nu_reactants'       , data=fillArray( reNuReactants, 0.0 ) )
        group.create_dataset( r'ord_reactants'      , data=fillArray( reOrdReactants, 0.0 ) )

        group.create_dataset( r'n_products'         , data=reNProducts )
        group.create_dataset( r'products'           , data=fillArray( reProducts, -1 ) )
        group.create_dataset( r'nu_products'        , data=fillArray( reNuProducts, 0.0 ) )
        group.create_dataset( r'ord_products'       , data=fillArray( reOrdProducts, 0.0 ) )

        group.create_dataset( r'sum_nu'             , data=reSumNu )

        group.create_dataset( r'has_rev_arr'        , data=reHasRevArr )
        group.create_dataset( r'rev_arr_coeff'      , data=fillArray( reRevArrCoeff, 0.0, maxLength=(lambda x: 3) ) )
        group.create_dataset( r'has_high_arr'       , data=reHasHighArr )
        group.create_dataset( r'has_low_arr'        , data=reHasLowArr )
        group.create_dataset( r'adv_arr_coeff'      , data=fillArray( reAdvArrCoeff, 0.0, maxLength=(lambda x: 3) ) )
        group.create_dataset( r'has_troe'           , data=reHasTroe )
        group.create_dataset( r'n_troe_coeff'       , data=reNTroeCoeff )
        group.create_dataset( r'troe_coeff'         , data=fillArray( reTroeCoeff, 0.0, maxLength=(lambda x: 4) ) )

        group.create_dataset( r'has_efficiencies'   , data=reHasEfficiencies )
        group.create_dataset( r'n_efficiencies'     , data=reNEfficiencies )
        group.create_dataset( r'sp_efficiencies'    , data=fillArray( reSpEfficiencies, -1 ) )
        group.create_dataset( r'efficiencies'       , data=fillArray( reEfficiencies, 0.0 ) )
