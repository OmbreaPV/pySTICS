import numpy as np

def harvest_criteria(teaugrain_list, mat, ind_debdes, coderecolte, codeaumin, h2ograinmin, h2ograinmax):
    ''''
    This module checks the harvest criteria chosen by user : physiological maturity or grain/fruit water content.
    '''
    index = np.array([i for i in range(len(teaugrain_list))])
    maturity_criteria_list = (coderecolte == 1) & (mat == 1) # maturité physiologique atteinte
    water_content_min_criteria_list = (coderecolte == 2) & (codeaumin == 1) & (index >= ind_debdes) & (teaugrain_list > h2ograinmin)
    water_content_max_criteria_list = (coderecolte == 2) & (codeaumin == 2) & (index >= ind_debdes) & (teaugrain_list < h2ograinmax)

    return maturity_criteria_list, water_content_min_criteria_list, water_content_max_criteria_list