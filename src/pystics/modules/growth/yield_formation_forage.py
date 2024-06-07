


def yield_formation_forage(masec, msresjaune, msneojaune, msresiduel):
    '''
    This module computes the harvestable pool (mafruit) and actually harvested pool (msrec_fou) for forage crops.
    The harvestable pool is the total aboveground biomass without the dead fallen pool (masec), from which we remove the newly produced senescent biomass (msneojaune) and the residual senescent biomass (msresjaune = mafeuiljaune).
    The actually harvested pool is the harvestable pool (mafruit) from which we remove the residual pool (msresiduel) after cutting.
    '''

    # Harvestable pool
    mafruit = masec - msresjaune - msneojaune

    # Harvested pool
    msrec_fou = mafruit - msresiduel

    return mafruit, msrec_fou
