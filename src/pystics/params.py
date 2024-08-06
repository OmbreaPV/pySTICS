from dataclasses import dataclass, field
import pandas as pd
import numpy as np
import xmltodict
import inspect
from unidecode import unidecode
import pkg_resources

PARAMS_FOLDER = pkg_resources.resource_filename(__name__, 'parametrization_files')

@dataclass
class CropParams:
    """
    This class contains all crop parameters.
    Default values can not be used as it is, but instantiation initializes all attributes with values from given or provided XML plant file.

    Initialization instructions :
        - To retrieve parameters from a custom plant XML file, just fill in its path in file_path attribute when the object is instantiated
        - If file_path is set to empty string (''), XML file stored in the xml_folder_path folder will be read based on the species and variety of the instance. If xml_folder_path is an empty string (''), the parametrization_files folder will be read, which contains files publicly realeased by STICS project team. 
        - Species and variety attributes must always been given when the object is instantiated.
        
    Warning with ZLABOUR / ZPENTE / ZPRLIM : they are not used in default XML files provided by STICS project team, so if these files are used, they need to be provided during instantiation.
    """

    # Mandatory attributes
    species: str = ''
    variety: str = ''

    # Paths to plant XML files
    file_path: str = ''
    xml_folder_path: str = ''

    # Parameters whose default values are important
    ZLABOUR: int = 15
    ZPENTE: int = 100
    ZPRLIM: int = 150

    CODEPLANTE: str = ''
    CODEPERENNE: int = 0
    CODEINDETERMIN: int = 0
    HERBACEOUS : bool = True
    CODEBFROID: int = 0
    CODEPHOT: int = 0
    CODCALINFLO: int  = 0
    CODEIR: int = 0
    CODETEMPRAC: int = 0
    CODEGERMIN: int = 0
    CODEHYPO: int = 0
    CODERETFLO: int = 0
    CODGELLEV: int = 0
    CODGELJUV: int = 0
    CODGELFLO: int = 0
    CODGELVEG: int = 0
    CODETEMP: int = 0
    CODETREMP: int = 0
    TCMIN: int = 0
    TCMAX: int = 0
    JVC: int = 0
    TFROID: float = 0.
    AMPFROID: int = 0
    JVCMINI: int = 0
    PHOBASE: float = 0.
    PHOSAT: int = 0
    SENSIPHOT: int = 0
    TDMINDEB: float = 0.
    TDMAXDEB: float = 0.
    Q10: float = 0.
    STDORDEBOUR: int = 0
    STPLTGER: int = 0
    STLEVAMF: int = 0
    STAMFLAX: int = 0
    STDRPMAT: int = 0
    STLEVDRP: int = 0
    STFLODRP: int = 0
    STDRPNOU: int = 0
    STDRPDES: int = 0
    STLAXSEN: int = 0
    STSENLAN: int = 0
    CROIRAC: float = 0.
    KMAX: float = 0.
    ADENS: float = 0.
    BDENS: int = 0
    DLAIMAXBRUT: float = 0.
    DLAIMAX: float = 0.
    CODLAINET: int = 0
    PENTLAIMAX: float = 0.
    VLAIMAX: float = 0.
    UDLAIMAX: float = 0.
    LAICOMP: float = 0.
    DLAIMIN: float = 0.
    EXTIN: float = 0.
    TEOPT: float = 0.
    TEOPTBIS: float = 0.
    TEMIN: float = 0.
    TEMAX: float = 0.
    EFCROIJUV: float = 0.
    EFCROIVEG: float = 0.
    EFCROIREPRO: float = 0.
    PSISTO: float = 0.
    RAYON: float = 0.
    TDEBGEL: float = 0.
    TLETALE: float = 0.
    TGELFLO10: float = 0.
    TGELFLO90: float = 0.
    TGELJUV10: float = 0.
    TGELJUV90: float = 0.
    TGELVEG10: float = 0.
    TGELVEG90: float = 0.
    TGELLEV10: float = 0.
    TGELLEV90: float = 0.
    NBFGELLEV: int = 0
    PHYLLOTHERME: int = 0
    PGRAINMAXI: float = 0.
    NBJGRAIN: int = 0
    CGRAINV0: float = 0.
    CGRAIN: float = 0.
    NBGRMAX: int = 0
    NBGRMIN: int = 0
    VITIRCARB: float = 0.
    VITIRCARBT: float = 0.
    IRMAX: float = 0.
    TMINREMP: float = 0.
    TMAXREMP: float = 0.
    AFRUITPOT: float = 0.
    NBINFLO: int = 0
    INFLOMAX: int = 0
    PENTINFLORES: float = 0.
    SPFRMIN: float = 0.
    SPFRMAX: float = 0.
    DUREEFRUIT: int = 0
    NBOITE: int = 0.
    ALLOCFRMAX: float = 0.
    CFPF: float = 0.
    DFPF: float = 0.
    AFPF: float = 0.
    BFPF: float = 0.
    TDMIN: float = 0.
    TDMAX: float = 0.
    TCXSTOP: float = 0.
    ALPHACO2: float = 0.
    SLAMAX: float = 0.
    TIGEFEUIL: float = 0.
    PROPRES: float = 0.
    PROPRESP: float = 0.
    SLAMIN: float = 0.
    ABSCISSION: float = 0.
    RATIOSEN : float = 0.
    DURVIEF: float = 0.
    RATIODURVIEI: float = 0.
    RAPSENTURG: float = 0.
    SWFACMIN: float = 0.
    STRESSDEV: float = 0.
    REMOBRES: float = 0.
    EFREMOBIL: float = 0.
    RESPLMAX: float = 0.
    PSITURG: float = 0.
    RSMIN: float = 0.
    STEMFLOWMAX: float = 0.
    MOUILLABIL: float = 0.
    KSTEMFLOW: float = 0.
    HAUTMAX: float = 0.
    HAUTBASE: float = 0.
    CODEBESO: int = 0
    CODEDORMANCE: int = 0
    CODEINTERCEPT: int = 0
    POTGERMI: float = 0.
    PROPJGERMIN: float = 0.
    TGMIN: float = 0.
    NBJGERLIM: int = 0
    SENSRSEC: float = 0.
    H2OFRVERT: float = 0.
    DESHYDBASE: float = 0.
    TEMPDESHYD: float = 0.
    ELMAX: float = 0.
    BELONG: float = 0.
    CELONG: float = 0.
    NLEVLIM1: int = 0.
    NLEVLIM2: int = 0.
    STOPRAC: str = ''
    CONTRDAMAX: float = 0.
    COEFMSHAUT: float = 0.
    KHAUT: float = 0.
    JULVERNAL: float = 0.
    ENVFRUIT: float = 0
    CODEPHOT_PART: int = 0
    TUSTRESSMIN: float = 0.
    DURVIESUPMAX: float = 0.
    STOPRAC: str = ''

    STOPFEUILLE: float = field(init=False)

    def __post_init__(self):

        # We store ZPRLIM / ZPENTE / ZLABOUR values given by user because they will be erased by values from XML
        zprlim = self.ZPRLIM
        zpente = self.ZPENTE
        zlabour = self.ZLABOUR
        
        # Retrieve XML folder
        if self.xml_folder_path == '':
            self.xml_folder_path = PARAMS_FOLDER + '/plant'
        
        if self.file_path == '':
            species_to_plant_files = {'wheat' : 'wheat_plt' if unidecode(self.variety.lower()) in ['arminda','talent','thesee','soissons','promentin','sideral','thetalent','thesarmin','shango'] else f'DurumWheat_{self.variety.upper()}_plt',
                                      'rapeseed' : 'rapeseed_plt',
                                      'corn' : 'corn_plt',
                                      'alfalfa' : 'proto_alfalfa_plt',
                                      'fescue' : 'proto_fescue_plt',
                                      'timothy' : 'timothy_plt',
                                      'barley':'proto_barley_plt',
                                      'sorghum':'proto_sorghum_plt',
                                      'soybean':'proto_soybean_plt',
                                      'sunflower':'proto_sunflower_plt',
                                      'winter_barley':'proto_winterbarley_plt',
                                      'pea':'pea_plt',
                                      'oat':'BristleOat_CoverCrop_plt',
                                      'clover':'CrimsonClover_CoverCrop_plt',
                                      'grass':'grass_plt',
                                      'ryegrass':'ryegrass_CoverCrop_plt',
                                      'mustard':'mustard_CoverCrop_plt',
                                      'alfalfa':'proto_alfalfa_plt',
                                      'fescue':'proto_fescue_plt',
                                      'flax':'proto_flax_plt',
                                      'timothy':'timothy_plt',
                                      'vesce':'vetch_CoverCrop_plt', # vesce (fève, féverolle)
                                        }
            self.file_path = self.xml_folder_path + f"/{species_to_plant_files[f'{self.species.lower()}']}.xml"

        # Retrieve parameters values from XML file.
        self.get_params_from_xml()

        # We retrieve ZPRLIM / ZPENTE / ZLABOUR values given in instantiation if they are not provided in XML file.
        if self.ZPRLIM == -999:
            self.ZPRLIM = zprlim
            self.ZPENTE = zpente
            self.ZLABOUR = zlabour

        # We modify CODEBFROID if JVC and JVCMINI values are not consistent
        if (self.CODEBFROID == 2) & (self.JVC < self.JVCMINI):
            self.CODEBFROID = 1
        
        # Computation of stopfeuille according to STICS code.
        if self.CODEINDETERMIN == 1:
            self.STOPFEUILLE = 'LAX'
        elif self.CODEINDETERMIN == 2:
            self.STOPFEUILLE = 'SEN'
        
        # Custom attribute for herbaceous / ligneous and forage crops
        self.HERBACEOUS = True if self.CODEPERENNE == 1 else False
        if (self.species.lower() in ['fescue', 'alfalfa', 'timothy', 'dactyle']) | (self.variety.lower() in ['dactyle']):
            self.HERBACEOUS = True
        
        if self.CODLAINET == 2:
            self.DLAIMAX = self.DLAIMAXBRUT

    def get_params_from_xml(self):
        ''' 
        This method reads plant file (xml) from file_path and stores parameters values in dataclass fields.
        '''
        
        # 0. XML parsing
        with open(self.file_path, 'r') as f:
            data = f.read()
        dico = xmltodict.parse(data)

        # 1. Extract parameters values
        results = list(set(gen_dict_extract('#text','@nom',dico))) # extract all params from xml
        attributes = [i[0] for i in self.__dict__.items()] # list of CropParams fields
        crop_attr = [i.lower() for i in attributes if '__' not in i]
        r = [i for i in results if i[0].lower() in crop_attr]
        r_attr = [i[0].upper() for i in r]

        # stoprac parameter
        self.__setattr__('STOPRAC',r[r_attr.index('STOPRAC')][1].upper())
        del r[r_attr.index('STOPRAC')]
        del r_attr[r_attr.index('STOPRAC')]

        # codeplante parameter
        self.__setattr__('CODEPLANTE',r[r_attr.index('CODEPLANTE')][1].upper())
        del r[r_attr.index('CODEPLANTE')]
        del r_attr[r_attr.index('CODEPLANTE')]
        
        # Set values in dataclass fields
        for i, attr in zip(r,r_attr):
            self.__setattr__(attr, num(i[1]))

        # 2. Extract options values
        results2 = list(set(gen_dict_extract('@choix','@nomParam',dico)))
        r = [i for i in results2 if i[0].lower() in crop_attr]
        for i in r:
            attr = i[0].upper()
            self.__setattr__(attr,num(i[1]))

        # 3. Extract variety-specific parameters
        if any(i in self.file_path for i in ['wheat_plt','corn_plt','rapeseed_plt', 'timothy_plt', 'proto_barley']):
            for i in dico['fichierplt']['formalisme'][14]['tv']['variete']:
                if unidecode(i['@nom'].lower()) == unidecode(self.variety.lower()):
                    results3 = list(set(gen_dict_extract('#text','@nom',i)))
                    break
            attributes = [i[0] for i in self.__dict__.items()]
            crop_attr = [i.lower() for i in attributes if '__' not in i]
            r = [i for i in results3 if i[0].lower() in crop_attr] 
            for i in r:
                attr = i[0].upper()
                self.__setattr__(attr,num(i[1]))

@dataclass
class SoilParams:
    """
    This class contains all soil parameters.

    Initialization instructions :
        - all parameters values can be provided when the object is instantiated, or one of the 3 options of soil_source can be used to fill the parameters.
        - soil_source has 3 possible options : 'xml_file', 'db', 'ptf'.
            - xml_file : retrieve parameters from XML file stored in file_path. HMINF and HCCF must be provided in % (between 0 and 100) in XML file in this case.
            - db : used if a xarray dataset is provided in ds field, containing all parameters values, and latitude/longitude are used to select the location in the dataset.
            - ptf : used if only CLAY/SILT/SAND/SOC/DAF_1 are provided, it computes other parameters values with pedotransfer functions described in STICS book.
    """

    source: str = ''
    ds = None # dataset with soil parameters in different locations in France
    file_path: str = ''
    soil_name: str = ''
    latitude: float = 0.
    longitude: float = 0.
    CLAY: float = 0.
    SILT: float = 0.
    SAND: float = 0.
    SOC: float = 0.
    ARGI: float = 0.
    DAF_1: float = 0.
    DAF_2: float = 0.
    DAF_3: float = 0.
    DAF_4: float = 0.
    DAF_5: float = 0.
    EPC_1: int = 0
    EPC_2: int = 0
    EPC_3: int = 0
    EPC_4: int = 0
    EPC_5: int = 0

    # Parameters computed with pedotransfer functions
    Q0: float = 0.
    HMINF_1: float = 0.
    HMINF_2: float = 0.
    HMINF_3: float = 0.
    HMINF_4: float = 0.
    HMINF_5: float = 0.
    HCCF_1: float = 0.
    HCCF_2: float = 0.
    HCCF_3: float = 0.
    HCCF_4: float = 0.
    HCCF_5: float = 0.
    
    # Constant parameters for all soils
    ALBEDO: float = 0.2 # TODO : compute soil texture from clay/silt/sand and then use table 15.24 from STICS documentation. Until then : albedo = 0.2
    ZESX: int = 60
    CFES: float = 5.0
    Z0SOLNU: float = 0.01

    DEPTH: int = field(init=False)
    HCC: np.ndarray = field(init=False)
    HMIN: np.ndarray = field(init=False)
    DAF: np.ndarray = field(init=False)

    def __post_init__(self):

        if self.source == 'xml_file':
            self.get_params_from_xml()
             
            # We convert variables linked to water content to mm water .cm soil-1
            self.HMINF_1 = self.HMINF_1 / 100 * 10
            self.HMINF_2 = self.HMINF_2 / 100 * 10
            self.HMINF_3 = self.HMINF_3 / 100 * 10
            self.HMINF_4 = self.HMINF_4 / 100 * 10
            self.HMINF_5 = self.HMINF_5 / 100 * 10
            self.HCCF_1 = self.HCCF_1 / 100 * 10
            self.HCCF_2 = self.HCCF_2 / 100 * 10
            self.HCCF_3 = self.HCCF_3 / 100 * 10
            self.HCCF_4 = self.HCCF_4 / 100 * 10
            self.HCCF_5 = self.HCCF_5 / 100 * 10

            # Field capacity and wilting point must be fixed by bulk density.
            self.HMINF_1 = self.HMINF_1 * self.DAF_1
            self.HMINF_2 = self.HMINF_2 * self.DAF_2
            self.HMINF_3 = self.HMINF_3 * self.DAF_3
            self.HMINF_4 = self.HMINF_4 * self.DAF_4
            self.HMINF_5 = self.HMINF_5 * self.DAF_5
            self.HCCF_1 = self.HCCF_1 * self.DAF_1
            self.HCCF_2 = self.HCCF_2 * self.DAF_2
            self.HCCF_3 = self.HCCF_3 * self.DAF_3
            self.HCCF_4 = self.HCCF_4 * self.DAF_4
            self.HCCF_5 = self.HCCF_5 * self.DAF_5

        elif self.source == 'db':
            self.get_params_from_netcdf()
        elif self.source == 'ptf':
            self.get_params_from_ptf()

        self.EPC_1 = int(self.EPC_1)
        self.EPC_2 = int(self.EPC_2)
        self.EPC_3 = int(self.EPC_3)
        self.EPC_4 = int(self.EPC_4)
        self.EPC_5 = int(self.EPC_5)
        self.DEPTH = int(self.EPC_1 + self.EPC_2 + self.EPC_3 + self.EPC_4 + self.EPC_5)
        self.ZESX = min(self.ZESX, self.DEPTH)
        self.HCC = np.array([self.HCCF_1 for i in range(self.EPC_1)] + [self.HCCF_2 for i in range(self.EPC_2)] + [self.HCCF_3 for i in range(self.EPC_3)] + [self.HCCF_4 for i in range(self.EPC_4)] + [self.HCCF_5 for i in range(self.EPC_5)])
        self.HMIN = np.array([self.HMINF_1 for i in range(self.EPC_1)] + [self.HMINF_2 for i in range(self.EPC_2)] + [self.HMINF_3 for i in range(self.EPC_3)] + [self.HMINF_4 for i in range(self.EPC_4)] + [self.HMINF_5 for i in range(self.EPC_5)])
        self.DAF = np.array([self.DAF_1 for i in range(self.EPC_1)] + [self.DAF_2 for i in range(self.EPC_2)] + [self.DAF_3 for i in range(self.EPC_3)] + [self.DAF_4 for i in range(self.EPC_4)] + [self.DAF_5 for i in range(self.EPC_5)])

    def get_params_from_ptf(self):
        '''
        Pedotransfert functions are used to compute q0 / hmin / hccf
        '''
        self.Q0 = ptf_q0(self.CLAY, self.SILT, self.SAND)
        self.HMINF_1 = ptf_hminf(self.CLAY, self.SAND, self.SOC, self.DAF_1)
        self.HCCF_1 = ptf_hccf(self.CLAY, self.SAND, self.SOC, self.DAF_1)
        self.ARGI = self.CLAY
    
    def get_params_from_xml(self):
        '''
        This function retrieves soil parameters values from XML file.
        '''

        with open(self.file_path, 'r') as f:
            data = f.read()
        dico = xmltodict.parse(data)
        for dic in dico['sols']['sol']:
            if dic['@nom'] == self.soil_name:
                results = list(set(gen_dict_extract('#text','@nom',dic)))
                for layer in dic['tableau']:

                    r = list(set(gen_dict_extract('#text','@nom',layer)))
                    for attr in ['DAF_','HMINF_','HCCF_','EPC_']:
                        attribut = attr+f"{layer['@nom'][-1]}"
                        if attr == 'EPC_':
                            self.__setattr__(attribut, float([i[1] for i in r if i[0] == attr[:-1].lower()][0]))
                        else:
                            self.__setattr__(attribut, float([i[1] for i in r if i[0] == attr[:-1]][0]))
        
        attributes = [i[0] for i in self.__dict__.items()]
        soil_attr = [i.lower() for i in attributes if '__' not in i]

        r = [i for i in results if i[0].lower() in soil_attr]

        for i in r:
            attr = i[0].upper()
            self.__setattr__(attr,num(i[1]))

        self.ZESX = int(self.ZESX)

    def get_params_from_netcdf(self):
        '''
        This function retrieves soil parameters values from a xarray dataset based on location. The soil considered has only one layer (parameters are homogeneous on whole depth).
        Xarray is necessary to run this method (pip install xarray)
        '''

        selec = self.ds.sel(x=self.longitude, y=self.latitude, method = "nearest")

        for attr in ['EPC','DAF','HMINF','HCCF']:
            for layer in range(1,6):
                attribut_name = attr+"_"+str(layer)

                if layer == 1:
                    self.__setattr__(attribut_name, float(selec[attr].values))
                else :
                    self.__setattr__(attribut_name, 0)

        self.Q0 = float(selec["Q0"].values)
        self.ARGI = float(selec["clay"].values)
        self.SAND = float(selec["sand"].values)
        self.SILT = float(selec["silt"].values)
        self.CLAY = float(selec["clay"].values)


@dataclass
class Constants:
    """
    This class contains all constants parameters.
    """

    COEFB: float = (0.0815 / 100)  # For homogenity with ebmax (in cg.MJ-1), we divide coefb (in g.MJ-1) by 100
    PROPRAC: float = 0.2
    Y0MSRAC: float = 0.7
    LVOPT: float = 0.5
    PARSURRG: float = 0.48
    DIFTHERM: float = 0.00537
    PSIHUCC: float = -0.03
    PSIHUMIN: float = -1.5
    IRRLEV: float = 20.
    DASEUILBAS: float = 1.4
    DASEUILHAUT: float = 2.


@dataclass
class ManageParams:
    """
    This class contains crop management parameters.

    Initialization instructions : 
        - To use a XML file to instantiate ManageParams, a file_path must be provided.
    """
    
    file_path: str = ''
    IPLT0: int = 0
    PROFSEM: float = 0
    DENSITESEM: int = 0
    CODRECOLTE: int = 0
    CODEAUMIN: int = 0
    H2OGRAINMAX: float = 0.
    H2OGRAINMIN: float = 0.
    CODLOCIRRIG: int = 0
    CODECALIRRIG: int = 0
    CODEDATE_IRRIGAUTO: int = 0
    DATEDEB_IRRIGAUTO: int = 0
    DATEFIN_IRRIGAUTO: int = 0
    DOSIMX: int = 0
    DOSEIRRIGMIN: int = 0
    RATIOL: float = 0.
    EFFIRR: float = 0.
    CODEDATEAPPH2O: int = 0
    IRRIGATION_INTERVENTIONS: dict = field(default_factory= lambda: {})
    TEMPFAUCHE: float = 0.
    CODEMODFAUCHE: int = 0
    HAUTCOUPE: float = 0.
    CODEFAUCHE: int = 0
    JULFAUCHE: list = field(init=False)
    MSCOUPEMINI: list = field(init=False)
    

    def get_params_from_xml(self):
        
        with open(self.file_path, 'r') as f:
            data = f.read()
        dico = xmltodict.parse(data)
        results = list(set(gen_dict_extract('#text','@nom',dico)))
        list_attr = [i[0] for i in results]
        
        for i in inspect.getmembers(self):
            attr = i[0]
            if attr.lower() in list_attr:
                self.__setattr__(attr, float([i[1] for i in results if i[0] == attr.lower()][0]))

        results2 = list(set(gen_dict_extract('@choix','@nomParam',dico)))
        list_attr = [i[0].lower() for i in results2]

        for i in inspect.getmembers(self):
            attr = i[0]
            if attr.lower() in list_attr:
                self.__setattr__(attr, float([i[1] for i in results2 if i[0].lower() == attr.lower()][0]))

        if (self.CODECALIRRIG == 2) & (int(dico['fichiertec']['formalisme'][4]['option'][0]['choix'][1]['ta']['@nb_interventions']) > 0): # si nb_interv > 0
            if int(dico['fichiertec']['formalisme'][4]['option'][0]['choix'][1]['ta']['@nb_interventions']) == 1:
                self.IRRIGATION_INTERVENTIONS[f"{dico['fichiertec']['formalisme'][4]['option'][0]['choix'][1]['ta']['intervention']['colonne'][0]['#text']}"] = float(dico['fichiertec']['formalisme'][4]['option'][0]['choix'][1]['ta']['intervention']['colonne'][1]['#text'])
            else:
                for i in dico['fichiertec']['formalisme'][4]['option'][0]['choix'][1]['ta']['intervention']:
                    self.IRRIGATION_INTERVENTIONS[f"{i['colonne'][0]['#text']}"] = float(i['colonne'][1]['#text'])

        if (self.CODEFAUCHE == 1):
            if (int(dico['fichiertec']['formalisme'][7]['option'][0]['choix'][0]['option'][2]['choix'][1]['ta']['@nb_interventions']) > 0):
                if int(dico['fichiertec']['formalisme'][7]['option'][0]['choix'][0]['option'][2]['choix'][1]['ta']['@nb_interventions']) == 1:
                    self.JULFAUCHE = [int(dico['fichiertec']['formalisme'][7]['option'][0]['choix'][0]['option'][2]['choix'][1]['ta']['intervention']['colonne'][0]['#text'])]
                    self.MSCOUPEMINI = [int(dico['fichiertec']['formalisme'][7]['option'][0]['choix'][0]['option'][2]['choix'][1]['ta']['intervention']['colonne'][8]['#text'])]
                else:
                    l_julfauche = []
                    l_mscoupemini = []
                    for i in dico['fichiertec']['formalisme'][7]['option'][0]['choix'][0]['option'][2]['choix'][1]['ta']['intervention']:
                        l_julfauche.append(int(i['colonne'][0]['#text']))
                        l_mscoupemini.append(int(i['colonne'][8]['#text']))
                    self.JULFAUCHE = l_julfauche.copy()
                    self.MSCOUPEMINI = l_mscoupemini.copy()
        else:
            self.JULFAUCHE = []
            self.MSCOUPEMINI = []
    
    def __post_init__(self):
        if self.file_path != '':
            self.get_params_from_xml()

@dataclass
class StationParams:
    '''
    This class contains station parameters.

    Initialization instructions : 
        - To use a XML file to instantiate StationParams, a file_path must be provided.
    '''
    
    file_path: str = ''    
    CODEETP: int = 0
    LATITUDE: float = 0.
    CODERNET: int = 0
    CODECALTEMP: int = 0

    # constant parameters
    ALBVEG: float = 0.23
    ACLIM: float = 20
    AANGST: float = 0.18
    BANGST: float = 0.62
    CORECTROSEE: float = 1.
    ZR: float = 2.5
    ALPHAPT : float = 1.26

    def get_params_from_xml(self):
        with open(self.file_path, 'r') as f:
            data = f.read()
        dico = xmltodict.parse(data)
        results = list(set(gen_dict_extract('#text','@nom',dico)))

        attributes = [i[0] for i in self.__dict__.items()]
        station_attr = [i.lower() for i in attributes if '__' not in i]
        r = [i for i in results if i[0].lower() in station_attr]
        for i in r:
            attr = i[0].upper()
            self.__setattr__(attr,num(i[1]))

        results2 = list(set(gen_dict_extract('@choix','@nomParam',dico)))
        self.CODEETP = int([i[1] for i in results2 if i[0] if i[0] == 'codeetp'][0])
        self.CODERNET = int([i[1] for i in results2 if i[0] if i[0] == 'codernet'][0])
        self.CODECALTEMP = int([i[1] for i in results2 if i[0] if i[0] == 'codecaltemp'][0])
    
    def __post_init__(self):
        if self.file_path != '':
            self.get_params_from_xml()


@dataclass
class InitialParams:
    '''
    This class contains initial values.

    Initialization instructions : 
        - To use a XML file to instantiate InitialParams, a file_path must be provided.
        - If HINITF values unit is % between 0 and 100 (instead of mm water . cm soil-1), fill hinitf_unit = '%'.
    '''

    file_path: str = ''
    hinitf_unit: str = ''
    HINITF_1: float = 0. 
    HINITF_2: float = 0. 
    HINITF_3: float = 0. 
    HINITF_4: float = 0. 
    HINITF_5: float = 0.25
    DENSINITIAL_1: float = 0. 
    DENSINITIAL_2: float = 0. 
    DENSINITIAL_3: float = 0. 
    DENSINITIAL_4: float = 0. 
    DENSINITIAL_5: float = 0.25
    STADE0 : str = ''
    LAI0 : float = 0.
    ZRAC0 : float = 0.
    MASEC0 : float = 0.
    RESTEMP0: float = 0.

    def get_params_from_xml(self):
        with open(self.file_path, 'r') as f:
            data = f.read()
        dico = xmltodict.parse(data)
        for i in range(5):
            attr = f'HINITF_{i+1}'
            self.__setattr__(attr, float(dico['initialisations']['sol']['Hinitf']['horizon'][i]['#text']))
            attr = f'DENSINITIAL_{i+1}'
            self.__setattr__(attr, float(dico['initialisations']['plante'][0]['densinitial']['horizon'][i]['#text']))
        self.ZRAC0 = float(dico['initialisations']['plante'][0]['zrac0'])
        self.LAI0 = float(dico['initialisations']['plante'][0]['lai0'])
        self.STADE0 = dico['initialisations']['plante'][0]['stade0']
        num = int(dico['initialisations']['plante'][0]['option']['@choix'])
        self.MASEC0 = float(dico['initialisations']['plante'][0]['option']['choix'][num-1]['masec0'])
        
              
    def __post_init__(self):
        if self.file_path != '':
            self.get_params_from_xml()

        if self.hinitf_unit == '%' :
            # We convert variables linked to water content to mm water .cm soil-1 (given in % between 1 and 100)
            self.HINITF_1 = self.HINITF_1 / 100 * 10
            self.HINITF_2 = self.HINITF_2 / 100 * 10
            self.HINITF_3 = self.HINITF_3 / 100 * 10
            self.HINITF_4 = self.HINITF_4 / 100 * 10
            self.HINITF_5 = self.HINITF_5 / 100 * 10

def parametrization_from_stics_example_files(species, variety, xml_folder_path = PARAMS_FOLDER):
    '''
    This function allows to read parameters files from USM examples given by STICS project team.
    It reads weather/plant/soil/manage/station/initial files for the USM associated to species and variety (association is based on usms.xml file).
    It returns dataclass instances containing the parameters for each file.
    If not xml_folder_path is provided, the considered folder is parametrization_files (containing STICS project team example files).
    '''

    # Retrieve USM associated to species and variety chosen by user
    with open(xml_folder_path + '/example/usms.xml', 'r') as f:
        data = f.read()
    dico = xmltodict.parse(data)
    for dic in dico['usms']['usm']:
        if dic['@nom'].lower() == (species.lower() if species.lower() != 'corn' else 'maize'):
            usm_files = dic

    # Weather file
    weather = weather_parametrization(xml_folder_path=xml_folder_path, clim1=usm_files['fclim1'], clim2=usm_files['fclim2'])

    # Plant parameters
    crop = CropParams(species = species, variety = variety, xml_folder_path = xml_folder_path+'/plant')
    
    # Crop management parameters
    manage = ManageParams(file_path = xml_folder_path + f"/example/{usm_files['plante'][0]['ftec']}")

    # Soil parameters
    soil = SoilParams(file_path = xml_folder_path + '/example/sols.xml',
                      soil_name = f"{usm_files['nomsol']}",
                      source = 'xml_file')

    # Station parameters
    station = StationParams(file_path = xml_folder_path + f"/example/{usm_files['fstation']}")

    # Constants parameters
    constants = Constants()

    # Initial parameters
    initial = InitialParams(file_path = xml_folder_path + f"/example/{usm_files['finit']}", hinitf_unit = '%')
    
    return weather, crop, manage, soil, station, constants, initial


def weather_parametrization(xml_folder_path, clim1='', clim2='', ):
    '''
    This function reads weather files provided by STICS project team and returns a pandas dataframe in the right format for run_pystics_simulation.
    '''

    weather_year1 = pd.read_csv(xml_folder_path + f"/example/{clim1}", header=None).rename(columns={0:"raw"})
    weather_year1["raw"] = weather_year1["raw"].str.split()
    weather_year1 = pd.DataFrame(weather_year1.raw.tolist(), index= weather_year1.index, columns = ['file','year','month','day','doy','temp_min','temp_max','radiation','etp','rain','wind','tpm','co2'])

    if clim1 == clim2: # cas où 1 seul fichier climat = plante semée et récoltée sur la même année
        weather = weather_year1
    else:
        weather_year2 = pd.read_csv(xml_folder_path + f"/example/{clim2}", header=None).rename(columns={0:"raw"})
        weather_year2["raw"] = weather_year2["raw"].str.split()
        weather_year2 = pd.DataFrame(weather_year2.raw.tolist(), index= weather_year2.index, columns = ['file','year','month','day','doy','temp_min','temp_max','radiation','etp','rain','wind','tpm','co2'])

        weather = pd.concat([weather_year1, weather_year2], axis=0)

    weather = weather.reset_index()
    weather['date'] = pd.to_datetime(dict(year=weather.year, month=weather.month, day=weather.day))
    weather = weather.drop(['index','file'], axis=1)
    for col in ['doy','temp_min','temp_max','radiation','etp','rain','wind','tpm','co2']:
        weather[col] = weather[col].astype('float')

    return weather

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

def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)
    

def ptf_q0(clay, silt, sand):
    '''
    Pedotransfer function to compute q0.
    See M-R. Dobarco et al. 2019. Pedotransfer functions for predicting available water capacity in French soils, their applicability domain and associated uncertainty. Geoderma.
    '''
    return np.where((clay >= 50) & (silt <= 50),
                    5+0.06*(100-clay),
                    np.where((clay <= 20) & (sand >= 80),
                            5 + 0.15 * (100 - sand),
                            8 + 0.8 * clay))

def ptf_hminf(clay, sand, soc, daf):
    '''
    Pedotransfer function to compute hminf (mm water . cm soil-1).
    See M-R. Dobarco et al. 2019. Pedotransfer functions for predicting available water capacity in French soils, their applicability domain and associated uncertainty. Geoderma.
    '''
    return (-0.029 + 4.35 * 1e-3 * clay + -6.08 * 1e-5 * sand + 1.7 * 1e-2*  soc + 4.77 * 1e-2 * daf) * 10 # *10 : % to mm.cm-2

def ptf_hccf(clay, sand, soc, daf):
    '''
    Pedotransfer function to compute hccf (mm water . cm soil-1).
    See M-R. Dobarco et al. 2019. Pedotransfer functions for predicting available water capacity in French soils, their applicability domain and associated uncertainty. Geoderma.
    '''
    return (0.127 + 2.29 * 1e-3 * clay + -1.21 * 1e-3 * sand + 4.35 * 1e-2 * soc + 7.35 * 1e-2 * daf) * 10 # *10 : % to mm.cm-2
    