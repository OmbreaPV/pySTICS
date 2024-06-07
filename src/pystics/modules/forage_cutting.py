


def cutting(codefauche, cut_number, julfauche, doy, codemodfauche, msresiduel, masec, masecneo, lairesiduel, lai, msrec_fou, mafruit, mscoupemini, masectot, msresjaune, msneojaune, mafeuiljaune, msres, stade0, lev, amf, lax, flo, drp, debdes):
    '''

    This modules checks the 4 conditions for forage crops cutting, and triggers cutting if they are verified.
    
    Simplifications compared to STICS :
        - only codemodfauche = 2 implemented
        - only case where hautcoupe is an input is implemented
        - roots not affected by cutting
        - phenological stage after cutting is stade0
    '''
    
    if (codefauche == 1) & (cut_number < len(julfauche)):

        # When a cutting date has been reached, we compute the three other conditions before triggering cutting
        if (doy >= julfauche[cut_number]) & (codemodfauche == 2):
        
            # 3 conditions to trigger cutting
            cut_condition_lai = lai > lairesiduel
            cut_condition_masec = masec > msresiduel
            cut_condition_msrec_fou = msrec_fou > mscoupemini[cut_number]
            
            # If the 3 conditions are verified, we trigger the cutting
            if cut_condition_lai & cut_condition_masec & cut_condition_msrec_fou:

                cut_number += 1
                masectot = masectot + masec
                msrec_fou = 0
                mafruit = msresiduel
                masecneo = 0
                msresjaune = 0
                msneojaune = 0
                mafeuiljaune = 0 # --> mafeuiljaune becomes the senescent pool of newly produced biomass
                masec = msresiduel
                msres = msresiduel
                lai = lairesiduel
                lev = 1
                amf = 1 if stade0 == 'amf' else 0
                lax = 0
                flo = 0
                drp = 0
                debdes = 0            
    
    return masectot, msrec_fou, mafruit, masecneo, msresjaune, msneojaune, mafeuiljaune, masec, msres, lai, cut_number, lev, amf, lax, flo, drp, debdes

