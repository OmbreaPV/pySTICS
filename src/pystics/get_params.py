import xmltodict
import pandas as pd
from pystics.params import StationParams, CropParams, ManageParams, SoilParams, Constants, UserOptions
from pystics.meteo import get_cached_meteo

def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)


def parametrization_from_usm_example(species, variety, meteo_source, lat=0,lon=0,year=0):
    # '''
    # Input : un USM ou plutôt une espèce à laquelle on a associé en dur (dans dico) des fichiers ini/sol/climat/tec
    # Output : retourne les instances des dataclass utilisées en simulation avec les paramètres qui viennent des fichiers ini/sol/lcimat/tec
    # '''

    # Dico associant 1 espèce/variété aux fichiers plant/sol/climat/tec associés
    list_varieties_wheat_plt = ['Arminda','Talent','Thesee','Soissons','Promentin','Sideral','Thétalent','Thésarmin','Shango'] # liste des variétés dans wheat_plt.xml
    wheat_plt = False # sert à condition dans fonction param_from_plant_file
    if species == 'wheat':
        if variety in list_varieties_wheat_plt:
            wheat_plt = True # sert à condition dans fonction param_from_plant_file
            species_to_example_files = {'wheat':{'plant':f'plant_files/wheat_plt.xml', 'ini':'examples/ble_ini.xml', 'tec':'examples/Ble_tec.xml', 'soil':'solble', 'meteo_sta':'examples/climblej_sta.xml', 'meteo_year1':'examples/climblej.1994', 'meteo_year2':'examples/climblej.1995'}}
        else:
            species_to_example_files = {'wheat':{'plant':f'plant_files/DurumWheat_{variety.upper()}_plt.xml', 'ini':'examples/ble_ini.xml', 'tec':'examples/Ble_tec.xml', 'soil':'solble', 'meteo_sta':'examples/climblej_sta.xml', 'meteo_year1':'examples/climblej.1994', 'meteo_year2':'examples/climblej.1995'}}
    

    ##############################
    ### Lecture des fichiers meteo
    if meteo_source == 'stics':
        meteo_year1 = pd.read_csv(f"../parametrization_files/{species_to_example_files[species]['meteo_year1']}", header=None).rename(columns={0:"raw"})
        meteo_year1["raw"] = meteo_year1["raw"].str.split()
        meteo_year1 = pd.DataFrame(meteo_year1.raw.tolist(), index= meteo_year1.index, columns = ['file','year','month','day','doy','Temp_min','Temp_max','Radiation','ETP','Rain','Wind','TPM','CO2'])

        meteo_year2 = pd.read_csv(f"../parametrization_files/{species_to_example_files[species]['meteo_year2']}", header=None).rename(columns={0:"raw"})
        meteo_year2["raw"] = meteo_year2["raw"].str.split()
        meteo_year2 = pd.DataFrame(meteo_year2.raw.tolist(), index= meteo_year2.index, columns = ['file','year','month','day','doy','Temp_min','Temp_max','Radiation','ETP','Rain','Wind','TPM','CO2'])

        meteo = pd.concat([meteo_year1,meteo_year2], axis=0)
        meteo = meteo.reset_index()
        meteo['Date'] = pd.to_datetime(dict(year=meteo.year, month=meteo.month, day=meteo.day))
        meteo = meteo.drop(['index','file'], axis=1)
        for col in ['doy','Temp_min','Temp_max','Radiation','ETP','Rain','Wind','TPM','CO2']:
            meteo[col] = meteo[col].astype('float')
        meteo = meteo.rename(columns={'year':'ANNEE', 'doy':'DOY'})
    # Pas de temp moy...

    elif meteo_source == 'era5':
        meteo = get_cached_meteo(lat,lon,year,year+1)

    #########################
    # Lecture du fichier tec
    manage = ManageParams()
    with open(f"../parametrization_files/{species_to_example_files[species]['tec']}", 'r') as f:
        data = f.read()
    dico = xmltodict.parse(data)
    results = list(set(gen_dict_extract('#text','@nom',dico)))
    manage.IPLT0 = float([i[1] for i in results if i[0] == 'iplt0'][0])
    manage.H2OGRAINMAX = float([i[1] for i in results if i[0] == 'h2ograinmax'][0])

    results2 = list(set(gen_dict_extract('@choix','@nomParam',dico)))
    manage.CODRECOLTE = float([i[1] for i in results2 if i[0] if i[0] == 'codrecolte'][0])
    manage.CODEAUMIN = float([i[1] for i in results2 if i[0] if i[0] == 'codeaumin'][0])

    ###########################
    # Lecture du fichier plante
    # Va chercher param plante + options (codedormance, codeintercept, codebeso)
    with open(f"../parametrization_files/{species_to_example_files[species]['plant']}", 'r') as f:
        data = f.read()
    if meteo_source=='stics':
        crop = CropParams(species=species, CO2=float(meteo.CO2[0]))
    else:
        crop = CropParams(species=species)
    crop = param_from_plant_file(data,crop, wheat_plt, variety) # pas même fonction exactement que celle qu'on lancera hors exemple STICS

    # Cas du blé : pour zprlim/zpente/zlabour ils sont à -999 parce que l'option dans wheat_plt.xml est true density, sauf que cette option n'est pas implémentée dans pyStics.
    if species == 'wheat': # je fais ça car j'ai pas encore regardé les autres espèces si j'avais le mm prob
        crop.ZPRLIM = 160 # fichier baresoil_plt
        crop.ZPENTE = 100 # fichier baresoil_plt
        crop.ZLABOUR = 20 # fichier baresoil_plt

    #########################
    # Fichier sol
    soil = SoilParams()
    with open(f"../parametrization_files/examples/sols.xml", 'r') as f:
        data = f.read()
    dico = xmltodict.parse(data)
    for dic in dico['sols']['sol']:
        if dic['@nom'] == species_to_example_files[species]['soil']:
            results = list(set(gen_dict_extract('#text','@nom',dic)))
            for layer in dic['tableau']:
                if layer['@nom'] == 'layer 1':
                    results_layer1 = list(set(gen_dict_extract('#text','@nom',layer)))
                if layer['@nom'] == 'layer 2':
                    results_layer2 = list(set(gen_dict_extract('#text','@nom',layer)))
            break
    # Commun à tout le sol
    attributes = [i[0] for i in soil.__dict__.items()] # on extrait la liste des attributs de la classe Crop
    soil_attr = [i.lower() for i in attributes if '__' not in i] # on met les noms des attributs en minuscule

    r = [i for i in results if i[0].lower() in soil_attr] # on prend l'intersection entre les attr de la classe Soil et les param trouvés dans le XML

    for i in r: # pour les param de l'intersection, en attribue à chaque attr la valeur extraite du XML
        attr = i[0].upper()
        setattr(soil,attr,num(i[1]))

    # Par layer
    soil.DAF_1 = float([i[1] for i in results_layer1 if i[0] == 'DAF'][0])
    soil.HMINF_1 = float([i[1] for i in results_layer1 if i[0] == 'HMINF'][0]) / 100
    soil.HCCF_1 = float([i[1] for i in results_layer1 if i[0] == 'HCCF'][0]) / 100
    soil.EPC_1 = int(float([i[1] for i in results_layer1 if i[0] == 'epc'][0]))

    soil.HMINF_2 = float([i[1] for i in results_layer2 if i[0] == 'HMINF'][0]) / 100
    soil.HCCF_2 = float([i[1] for i in results_layer2 if i[0] == 'HCCF'][0]) / 100
    soil.EPC_2 = int(float([i[1] for i in results_layer2 if i[0] == 'epc'][0]))

    # Variable de sol pas dans STICS :  hur1_init, hur2_init (pas de hur init tout court)



    #####################
    # Fichier station météo
    station = StationParams()
    with open(f"../parametrization_files/{species_to_example_files[species]['meteo_sta']}", 'r') as f:
        data = f.read()
    dico = xmltodict.parse(data)
    results = list(set(gen_dict_extract('#text','@nom',dico)))

    attributes = [i[0] for i in station.__dict__.items()] # on extrait la liste des attributs de la classe Crop
    station_attr = [i.lower() for i in attributes if '__' not in i] # on met les noms des attributs en minuscule

    r = [i for i in results if i[0].lower() in station_attr] # on prend l'intersection entre les attr de la classe Soil et les param trouvés dans le XML

    for i in r: # pour les param de l'intersection, en attribue à chaque attr la valeur extraite du XML
        attr = i[0].upper()
        setattr(station,attr,num(i[1]))

    results2 = list(set(gen_dict_extract('@choix','@nomParam',dico)))
    station.CODEETP = int([i[1] for i in results2 if i[0] if i[0] == 'codeetp'][0])


    # Quand on génère la météo, faut modifier dans station codeetp et latitude, les autres on garde par défaut
    if meteo_source=='era5':
        station.CODEETP = 2
        station.LATITUDE = lat # ça jsp trop

    # J'instancie Constants et User ici pour pas encombrer le Notebook
    constants = Constants()
    user = UserOptions()


    return meteo, crop, soil, manage, station, constants, user


def gen_dict_extract(key0, key1, var):
    if hasattr(var,'items'):
        for k, v in var.items():
            for k1, v1 in var.items():
                if (k == key0) & (k1==key1):
                    yield (v1,v)
                if isinstance(v, dict):
                    for result in gen_dict_extract(key0, key1, v):
                        yield result
                elif isinstance(v, list):
                    for d in v:
                        for result in gen_dict_extract(key0, key1, d):
                            yield result
                            
def param_from_plant_file(data, crop, wheat_plt, variety):
    ''' 
    Cas particulier du blé si variété dans heat_plt
    '''

    dico = xmltodict.parse(data)

    # Liste des attributs à ne pas chercher dans XML --> paramètres calculés ou qui n'existent pas de base dans STICS
    attr_not_searched = ['CO2', # donné input depuis autre classe
                            'FCO2', 'FCO2S','S','DURVIEI', # 
                            'DENSITESEM','ZRAC0','STADES','species','CODEPERENNE','CODEINDETERMIN', # entrés en dur pr chaque espèce
                            # 'CODEBFROID','CODEPHOT', 'IPLT0','CODCALINFLO','CODEIR' # options cherchées autrement dans XML ou dans autre XML
                            ]

    # 1. Extraction du XML des valeurs des attributs de la classe Crop
    results = list(set(gen_dict_extract('#text','@nom',dico))) # on extrait tous les noms de param et les valeurs qui existent dans XML
    attributes = [i[0] for i in crop.__dict__.items() if i[0] not in attr_not_searched] # on extrait la liste des attributs de la classe Crop
    crop_attr = [i.lower() for i in attributes if '__' not in i] # on met les noms des attributs en minuscule

    r = [i for i in results if i[0].lower() in crop_attr] # on prend l'intersection entre les attr de la classe Crop et les param trouvés dans le XML

    for i in r: # pour les param de l'intersection, en attribue à chaque attr la valeur extraite du XML
        attr = i[0].upper()
        setattr(crop,attr,num(i[1]))


    # 2. Extraction du XML des options de modélisation = pas si différent des attributs, mais pas même chemin dans XML
    results2 = list(set(gen_dict_extract('@choix','@nomParam',dico)))
    r = [i for i in results2 if i[0].lower() in crop_attr] # on prend l'intersection entre les attr de la classe Crop et les param trouvés dans le XML
    for i in r: # pour les param de l'intersection, en attribue à chaque attr la valeur extraite du XML
        attr = i[0].upper()
        setattr(crop,attr,num(i[1]))
    # crop.CODEBESO = float([i[1] for i in results2 if i[0] == 'codebeso'][0])
    # crop.CODEDORMANCE = float([i[1] for i in results2 if i[0] if i[0] == 'codedormance'][0])
    # crop.CODEINTERCEPT = float([i[1] for i in results2 if i[0] if i[0] == 'codeintercept'][0])

    if wheat_plt:
        # dans ce cas il faut aller chercher une 2ème fois une partie des param : ceux qui sont spécifiques à la variété. Les codes sont communs, pas besoin.
        # + prob avec les codes : ils ne sont pas accessibles facilement donc faut mettre des conditions en fonction des valeurs retournées
        for i in dico['fichierplt']['formalisme'][14]['tv']['variete']:
            if i['@nom'] == variety:
                results3 = list(set(gen_dict_extract('#text','@nom',i)))
                break
        
        attributes = [i[0] for i in crop.__dict__.items() if i[0] not in attr_not_searched] # on extrait la liste des attributs de la classe Crop
        crop_attr = [i.lower() for i in attributes if '__' not in i] # on met les noms des attributs en minuscule

        r = [i for i in results3 if i[0].lower() in crop_attr] # on prend l'intersection entre les attr de la classe Crop et les param trouvés dans le XML

        for i in r: # pour les param de l'intersection, en attribue à chaque attr la valeur extraite du XML
            attr = i[0].upper()
            setattr(crop,attr,num(i[1]))

    

    ## CODEBFROID
    codebfroid = float([i[1] for i in results2 if i[0] == 'codebfroid'][0])
    if codebfroid == 2: # 1 = pas de cold requir, 2 = besoin annuelle, 3 = besoin pérennes (mais aucun param nécessaire dans XML)
        if crop.JVC > crop.JVCMINI:
            crop.CODEBFROID = True
        elif crop.JVC < crop.JVCMINI: # condition dans le code source de STICS
            crop.CODEBFROID = False
    elif codebfroid == 3:
        crop.CODEBFROID = True
    else:
        crop.CODEBFROID = False

    ## CODEPHOT
    codephot = float([i[1] for i in results2 if i[0] == 'codephot'][0])
    if codephot == 1:
        crop.CODEPHOT = True
    elif codephot == 2:
        crop.CODEPHOT = False

    ## CODCALINFLO
    if crop.CODEINDETERMIN == 2: # concerne que plantes indéter
        codcalinflo = float([i[1] for i in results2 if i[0] if i[0] == 'codcalinflo'][0])
        if codcalinflo == 1:
            crop.CODCALINFLO = 1
        elif codcalinflo == 2:
            crop.CODCALINFLO == 2


    ## CODERETFLO = effet de l'eau sur le retard phasique avant DRP
    coderetflo = float([i[1] for i in results2 if i[0] if i[0] == 'coderetflo'][0])
    if coderetflo == 2: # non = pas de retard phasique = plante insensible à ce stress --> on fixe STRESSDEV à 0
        crop.STRESSDEV = 0
    # Rq : si coderetflo = 1 alors stressdev n'est pas à -999.

    return  crop


def param_from_plant_file_with_check_lists(data,crop): 

    dico = xmltodict.parse(data)

    # Liste des attributs à ne pas chercher dans XML --> paramètres calculés ou qui n'existent pas de base dans STICS
    attr_not_searched = ['CO2', # donné input depuis autre classe
                            'FCO2', 'FCO2S','S','DURVIEI', # 
                            'DENSITESEM','ZRAC0','STADES','species','CODEPERENNE','CODEINDETERMIN', # entrés en dur pr chaque espèce
                            # 'CODEBFROID','CODEPHOT', 'IPLT0','CODCALINFLO','CODEIR' # options cherchées autrement dans XML ou dans autre XML
                            ]

    # Extraction du XML des valeurs des attributs de la classe Crop
    results = list(set(gen_dict_extract('#text','@nom',dico))) # on extrait tous les noms de param et les valeurs qui existent dans XML
    attributes = [i[0] for i in crop.__dict__.items() if i[0] not in attr_not_searched] # on extrait la liste des attributs de la classe Crop
    crop_attr = [i.lower() for i in attributes if '__' not in i] # on met les noms des attributs en minuscule

    r = [i for i in results if i[0].lower() in crop_attr] # on prend l'intersection entre les attr de la classe Crop et les param trouvés dans le XML

    not_found = [i.upper() for i in crop_attr if i not in [i[0].lower() for i in r]] # not_found = attr de la classe Crop qui ne sont pas dans l'intersection

    for i in r: # pour les param de l'intersection, en attribue à chaque attr la valeur extraite du XML
        attr = i[0].upper()
        setattr(crop,attr,num(i[1]))

    param_999 = []

    # Extraction du XML des options de modélisation = pas si différent des attributs, mais pas même chemin dans XML
    results2 = list(set(gen_dict_extract('@choix','@nomParam',dico)))
    crop.CODEBESO = float([i[1] for i in results2 if i[0] == 'codebeso'][0])
    crop.CODEDORMANCE = float([i[1] for i in results2 if i[0] if i[0] == 'codedormance'][0])
    crop.CODEINTERCEPT = float([i[1] for i in results2 if i[0] if i[0] == 'codeintercept'][0])

    ## CODEBFROID
    codebfroid = float([i[1] for i in results2 if i[0] if i[0] == 'codebfroid'][0])
    if codebfroid == 2: # 1 = pas de cold requir, 2 = besoin annuelle, 3 = besoin pérennes (mais aucun param nécessaire dans XML)
        param_999 = param_999 + ['tdmindeb','tdmaxdeb'] # param qui sont à -999
        if crop.JVC > crop.JVCMINI:
            crop.CODEBFROID = True
        elif crop.JVC < crop.JVCMINI: # condition dans le code source de STICS
            crop.CODEBFROID = False
    elif codebfroid == 3:
        crop.CODEBFROID = True
        param_999 = param_999 + ['jvcmini','jvc','ampfroid','tfroid'] # param qui sont à -999 
    else:
        crop.CODEBFROID = False
        param_999 = param_999 + ['jvcmini','jvc','ampfroid','tfroid','tdmaxdeb','tdmindeb'] # param qui sont à -999 

    ## CODEPHOT
    crop.CODEPHOT = float([i[1] for i in results2 if i[0] if i[0] == 'codephot'][0])
    if crop.CODEPHOT == 1:
        crop.CODEPHOT = True
    elif crop.CODEPHOT == 2:
        crop.CODEPHOT = False
        param_999 = param_999 + ['sensiphot','phosat','phobase'] # param qui sont à -999

    ## CODEBESO
    crop.CODEBESO = float([i[1] for i in results2 if i[0] if i[0] == 'codebeso'][0])
    if crop.CODEBESO == 1: # approche coef K
        param_999 = param_999 + ['rsmin'] # resist stomat non utilisée dans cette approche.
    elif crop.CODEBESO == 2: # approche resistive
        param_999 = param_999 + ['kmax'] # coef de culture non utilisé dans cette approche.

    ## CODCALINFLO
    if crop.CODEINDETERMIN == 2: # concerne que plantes indéter
        crop.CODCALINFLO = float([i[1] for i in results2 if i[0] if i[0] == 'codcalinflo'][0])
        if crop.CODCALINFLO == 1:
            param_999 = param_999 + ['inflomax','pentinflores']
        elif crop.CODCALINFLO == 2:
            param_999 = param_999 + ['nbinflo']

    ## CODETRANSRAD
    codetransrad = float([i[1] for i in results2 if i[0] if i[0] == 'codetransrad'][0])
    if codetransrad == 1: # Beer law pr calculer interception rayonnem
        pass
    elif codetransrad == 2: # radiative transfer
        param_999 = param_999 + ['extin']

    ## CODERETFLO = effet de l'eau sur le retard phasique avant DRP
    coderetflo = float([i[1] for i in results2 if i[0] if i[0] == 'coderetflo'][0])
    if coderetflo == 2: # non = pas de retard phasique = plante insensible à ce stress --> on fixe STRESSDEV à 0
        crop.STRESSDEV = 0
        param_999 = param_999 + ['stressdev']
    # Rq : si coderetflo = 1 alors stressdev n'est pas à -999.

    ## CODEIR
    crop.CODEIR = float([i[1] for i in results2 if i[0] if i[0] == 'codeir'][0])
    if crop.CODEIR == 1: # utilise vitircarb pr calculer ratio grain/total biomass
        param_999 = param_999 + ['vitircarbt']
    elif crop.CODEIR == 2: # utilise vitircarbT pr calculer ratio grain/total biomass
        param_999 = param_999 + ['vitircarb']



    # Param à -999 selon si plante ann/pér et  déter / indéter
    if crop.CODEINDETERMIN == 1:
        param_999 = param_999 + ['afruitpot','inflomax','nbinflo','pentinflores','spfrmin','spfrmax','dureefruit','nboite','allocfrmax','cfpf','dfpf','afpf','bfpf', # param spécifiques au rempliss des fruits des plantes indéter
                                'stdrpnou' # pas de gdd drp - nou pour plantes déter
                                ] 

    elif crop.CODEINDETERMIN == 2:
        param_999 = param_999 + ['pgrainmaxi','nbjgrain','cgrainv0','cgrain','nbgrmax','nbgrmin','vitircarb','vitircarbt','irmax','tminremp','tmaxremp', # param spécifiques au rempliss des grains des plantes déter
                                'stdrpmat' # maturité pour les plantes pérennes calculée autrement
                                ]

    if crop.CODEPERENNE == 1:
        param_999 = param_999 + ['q10', # pas de dorm break pour plantes ann
                                'stdordebour' # pas de gdd dorm break - debour pour plantes ann
                                ]
    elif crop.CODEPERENNE == 2:
        param_999 = param_999 + ['stpltger'] # pas de gdd semis-germ pr plantes pérennes

    non_defined = [i[0] for i in [i for i in crop.__dict__.items() if i[0] not in attr_not_searched] if num(i[1]) <-98] # liste des attributs de la classe Crop dont les valeurs sont à -999
    # On ajoute une condition parce qu'il est normal que certains param soient à -999 selon si pér/annuelle ou codebfroid / codephot
    non_defined = [i for i in non_defined if i.lower() not in param_999]

    return  crop, not_found, non_defined

def initialize_meteo_and_param(date_semis, dormance, intercept_eau, species, year, lat, lon):

    '''
    Va chercher/construit la météo
    Crée les instances de classes crop / soil / constant / climate  / user
    Va chercher les param de la classe crop dans les XML
    '''

    # METEO
    meteo_raw = get_cached_meteo(lat,lon,year,year+1)

    # Instances crop / soil / constant / user
    dico_trad = {'vine':'vigne','wheat':'Ble'} # attention aux maj selon nom fichier _tec
    species_fr = dico_trad[species]
    soil = SoilParams()
    # climate = ClimateParams()
    crop = CropParams(species=species)
    constants = Constants()
    user = UserOptions()
    manage = ManageParams()

    # Prise en compte des choix user
    if date_semis == 'xml':
        with open(f'../parametrization_files/examples/{species_fr}_tec.xml', 'r') as f:
            data = f.read()
        dico = xmltodict.parse(data)
        results = list(set(gen_dict_extract('#text','@nom',dico)))
        manage.IPLT0 = float([i[1] for i in results if i[0] == 'iplt0'][0])
    else:
        manage.IPLT0 = int(date_semis)

    if dormance == 'calcul':
        user.CODEDORMANCE = 3 # CODEDORMANCE = 1 = fixée (IFINDORM à choisir en-dessous), 3 = date calculée
    elif dormance == 'fixe':
        user.CODEDORMANCE = 1
        user.IFINDORM = 390 # date sortie dormance choisie si CODEDORMANCE = 1

    if intercept_eau == 'oui':
        user.CODEINTERCEPT = 1 # 1 = yes, 2 = no. Si yes, il faut KSTEMFLOW, STEMFLOWMAX, MOUILLABIL
    elif intercept_eau == 'non':
        user.CODEINTERCEPT = 2


    # XML paramètres
    if species == 'vine':
        with open('../parametrization_files/plant_files/vine_GRENAC_plt.xml', 'r') as f:
            data = f.read()
    elif species == 'wheat':
        with open('../parametrization_files/plant_files/DurumWheat_BIENSUR_plt.xml', 'r') as f:
        # with open('../JavaSTICS-1.5.1-STICS-10.0.0/plant/DurumWheat_ACALOU_plt.xml', 'r') as f:
            data = f.read()


    crop, not_found, not_defined = param_from_plant_file_with_check_lists(data, crop)

    print('paramètres non trouvés : ',not_found)
    print('paramètres à -999 :  ',not_defined)


    # Valeurs dans les fichiers .tec qu'on ne va pas chercher. Ca doit dépendre des pratiques aussi, faudra faire une classe à part.
    if species == 'wheat':
        user.CODRECOLTE = 2 # choix du critère de récolte. 1 = physio maturity, 2 = water content
        user.CODEAUMIN = 2 # si CODERECOLTE = 2, choix de max/min. 1 = min, 2 = max.
        user.H2OGRAINMAX = 0.14 # contenu eau min pr récolte si CODERECOLTE = 2 et CODEAUMIN = 2

        # Param à -999 pour le blé
        crop.ZPRLIM = 160 # fichier baresoil_plt
        crop.ZPENTE = 100 # fichier baresoil_plt
        crop.ZLABOUR = 20 # fichier baresoil_plt

    elif species == 'vine':
        user.CODRECOLTE = 2 # choix du critère de récolte. 1 = physio maturity, 2 = water content
        user.CODEAUMIN = 2 # si CODERECOLTE = 2, choix de max/min. 1 = min, 2 = max.
        user.H2OGRAINMAX = 0.77 # contenu eau min pr récolte si CODERECOLTE = 2 et CODEAUMIN = 2

        # Param à -999 pour vigne
        crop.ZPRLIM = 150
        crop.ZPENTE = 100
        crop.ZLABOUR = 20 # fichier baresoil_plt
    
    return crop, soil, constants, user, manage, meteo_raw