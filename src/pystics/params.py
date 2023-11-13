from dataclasses import dataclass, field
from typing import Literal
import numpy as np

@dataclass
class CropParams:
    """Crop Parameters."""

    species: Literal[
        "wheat",
        "corn",
        "vine",
        "tomato",
        "rapeseed",
        "pea",
        "strawberry",
        "potato",
        "sorghum",
        "beet",
        "lettuce",
    ]

    CO2: float = 415  # CO2 concentration
    CODEBFROID: bool = False  # Cold requirement option: 1 = no cold requirement, 2 = vern days, 3 = dormancy with dorm break and post-dorm until bud break
    CODEPHOT: bool = False  # Is the crop sensitive to photoperiod?
    CODCALINFLO: int = 1  # 1 = number of inflorescences read in XML, 2 = number of inflorescences calculated
    CODEIR: int = 0  # Calculation ratio of grain to total biomass: 1 = proportional to time, 2 = proportional to thermal time

    # Plant Phenology
    TCMIN: int = 10  # Base temperature for growth [°C]
    TCMAX: int = 37  # Optimal temperature for growth [°C]

    # Vernalization (only if CODEBFROID = True)
    JVC: int = 45  # Required days of vernalization [days]
    TFROID: float = 6.5  # Optimal temperature for vernalization [°C]
    AMPFROID: int = 10  # Temperature range for vernalization [°C]
    JVCMINI: int = -999  # Minimum number of vernalizing days

    PHOBASE: float = 6.3  # Minimum day length to unlock photoperiod [hours]
    PHOSAT: int = 20  # Maximum day length above which there is no response to photoperiod [hours]
    SENSIPHOT: float = 0.8  # Photoperiod sensitivity: 1 = insensitive, 0.8 = very sensitive

    TDMINDEB: float = 5  # Minimum thermal threshold for hourly calculation of phasic duration between dormancy and bud breaks
    TDMAXDEB: float = 25  # Maximal temperature for hourly calculation of phasic duration between dormancy and bud breaks

    Q10: float = 2.17  # Threshold sum of temperature used for dormancy break

    # Plant Phasic Development
    STDORDEBOUR: int = 9145.4  # Cumulative time between dormancy break and bud break
    STPLTGER: int = 50  # Sum of degree days between sowing and germination
    ST2EMERGE: int = 20  # Sum of degree days between germination and emergence required to grow by 1 cm per day
    STLEVAMF: int = 225  # Sum of degree days between BBCH10 (emergence) and BBCH30 (1st node)
    STAMFLAX: int = 275  # Sum of degree days between BBCH30 (1st node) and BBCH39 (flag-Leaf)
    STDRPMAT: int = 725  # Sum of degree days between BBCH69 (DRP) and BBCH89 (maturity)
    STLEVDRP: int = 800  # Sum of degree days between BBCH10 (emergence) and BBCH69/71 (DRP)
    STFLODRP: int = 10  # Sum of degree days between BBCH65 (flowering) and BBCH69/71 (DRP)
    STDRPNOU: int = 91  # Sum of degree days between BBCH69 (DRP) and BBCH?? (NOU)

    # Plant growth parameters
    CROIRAC: float = 0.12  # Growth rate of the root front [cm/day]

    # Crop coefficient approach
    KMAX: int = 1  # Maximum stomatal conductance coefficient for water requirement (0 to 1)

    # Soil-related parameters
    ZLABOUR: int = 15  # Depth of tillage
    ZPENTE: int = 100  # Depth at which root density is 50% of surface root density
    ZPRLIM: int = 150  # Maximum root depth (reference)

    # LAI parameters
    ADENS: float = -0.65  # Interplant competition
    BDENS: int = 7  # Plant density threshold
    DLAIMAXBRUT: float = 0.00035  # Maximum rate of net LAI increase [m²/plant/growing degree days]
    PENTLAIMAX: float = 5.5  # Slope of the logistic function for LAI growth
    VLAIMAX: float = 1.8  # Ultimate LAI at the inflection point of the DELTAI = f(ULAI) function

    # Photosynthesis-related parameters
    EXTIN: float = 0.5  # Coefficient of PAR (Photosynthetically Active Radiation) extinction in the canopy
    TEOPT: float = 20  # Optimal temperature for growth (1/2)
    TEOPTBIS: float = 25  # Optimal temperature for growth (2/2)
    TEMIN: float = 10  # Minimum temperature for development
    EFCROIJUV: float = 1.2  # Maximum RUE (Radiation Use Efficiency) during the juvenile phase
    EFCROIVEG: float = 1.04  # Maximum RUE during the vegetative phase
    EFCROIREPRO: float = 2.25  # Maximum RUE during the reproductive phase

    # Stomatal close driving transpiration
    TIGEFEUIL: float = 0.5  # Variety-specific proportion of the total mass (stem + leaves) allocated to the stem
    PSISTO: float = 15  # Critical leaf water potential for stomatal closure
    RAYON: float = 0.02  # Average root radius

    # Frost stress index
    TDEBGEL: float = -1.5  # Temperature at which frost stress action begins
    TLETALE: float = -20  # Temperature at which the plant dies due to frost
    TGELFLO10: float = -2  # Temperature causing 10% frost damage on fruits during flowering
    TGELFLO90: float = -5  # Temperature causing 90% frost damage on fruits during flowering
    TGELJUV10: float = -2  # Temperature causing 10% frost damage during the juvenile phase
    TGELJUV90: float = -5  # Temperature causing 90% frost damage during the juvenile phase
    TGELVEG10: float = -2  # Temperature causing 10% frost damage during the vegetative phase
    TGELVEG90: float = -5  # Temperature causing 90% frost damage during the vegetative phase

    # Grain/fruit filling for DETERMINATE plants
    PGRAINMAXI: float = 3.33  # Maximum grain mass (at zero moisture content)
    NBJGRAIN: int = 10  # Number of days used to calculate the number of viable grains
    CGRAINV0: float = 10  # Fraction of the maximum number of grains when growth rate is zero
    CGRAIN: float = 0.1  # Slope of the function relating grain number to growth rate
    NBGRMAX: int = 10  # Maximum number of fruits per unit area
    NBGRMIN: int = 5  # Minimum number of fruits per unit area
    VITIRCARB: float = 0.0107  # Slope of the function relating progressive carbon HI to time since IDRP stage
    VITIRCARBT: float = 0  # Slope of the function relating progressive carbon HI to thermal time since IDRP stage
    IRMAX: float = 0.55  # Maximum harvest index (HI)
    TMINREMP: float = 0  # Temperature below which grain filling stops
    TMAXREMP: float = 37  # Temperature above which grain filling stops

    # Grain/fruit filling for INDETERMINATE plants
    AFRUITPOT: int = 3.66  # Maximum number of fruits per inflorescence per degree day
    NBINFLO: int = 15  # Number of inflorescences per plant
    INFLOMAX: int = 20  # Maximum number of inflorescences per plant
    PENTINFLORES: float = 2  # Parameter for calculating the number of inflorescences (optional)
    SPFRMIN: float = 0.75  # Minimum source-sink ratio for calculating trophic stress on fruit number
    SPFRMAX: float = 1  # Maximum source-sink ratio for calculating trophic stress on fruit number
    DUREEFRUIT: int = 276  # Duration from fruit formation to maturity
    NBOITE: int = 10  # Number of boxes/segments for age classes of fruits used to calculate fruit filling
    ALLOCFRMAX: float = 1  # Maximum daily allocation to fruits (unitless)

    CFPF: float = 15  # Parameter for the first phase of fruit growth (cell division)
    DFPF: float = 0.2  # Parameter for the first phase of fruit growth (cell division)
    AFPF: float = 0.55  # Logistic function parameter defining sink strength: relative age at which maximum growth occurs
    BFPF: float = 18  # Logistic function parameter defining sink strength: maximum growth relative to maximum fruit mass

    # Phenology
    TDMIN: float = 10  # Minimum temperature for development
    TDMAX: float = 37  # Maximum temperature for development
    TCXSTOP: float = 100  # Threshold temperature where leaf growth stops (used for phenology)

    # CO2 effect on RUE (Photosynthetic Efficiency) - an environmental parameter, not a plant parameter
    ALPHACO2: float = 1.2  # Sensitivity of growth to CO2 concentration

    ### BIOMASS PARTITIONING ###
    SLAMAX: float = 300  # Specific Leaf Area (SLA) - the ratio of leaf area to leaf mass
    TIGEFEUIL: float = 0.5  # Variety-specific proportion of mass allocated to the stem
    PROPRES: float = 0.1  # Maximum percentage of reserves that can be mobilized from above-ground organs
    PROPRESP: float = 0.1  # Maximum percentage of reserves that can be mobilized from storage organs of perennial plants
    SLAMIN: float = 180  # Minimum SLA of green leaves (cm²/g)
    ABSCISSION: float = 0  # Fraction of senescent leaves that fall to the ground

    # Senescence (MERLOT values)
    RATIOSEN: float = 0.8  # Proportion of leaf biomass lost due to senescence when the life duration is complete
    DURVIEF: float = 400  # Duration of the last produced leaves' life
    RATIODURVIEI: float = 0.8  # Ratio of the life duration of early leaves to the life duration of the last produced leaves
    RAPSENTURG: float = 0.05  # Threshold of active water content to simulate turgor stress senescence as a proportion of turgor stress

    # Water stress indices
    SWFACMIN: float = 0.1  # Minimum value of the water stress index
    STRESSDEV: float = 0.2  # Maximum phasic delay allowed by stresses

    # Reserves
    REMOBRES: float = 0.073  # Percentage of remobilizable C reserves per day
    EFREMOBIL: float = 0.4  # Efficiency of carbohydrate reserves remobilization
    RESTEMP0: float = 0  # Initial biomass of metabolic reserves in perennial organs
    RESPLMAX: float = 0.66  # Maximum biomass of temporary reserves

    INNSEN: float = 0.87  # Innate senescence rate

    # Threshold water content for reduced leaf growth
    PSITURG: float = 6  # Threshold water content limiting leaf growth (senescence for turfac calculation)

    # Fruits
    # ENVFRUIT: float = 0  # Proportion of pericarp among the maximum fruit mass

    # Energy budget
    RSMIN: float = 250  # Minimum stomatal resistance of leaves

    # Water intercepted by leaves
    STEMFLOWMAX: float = 0.3  # Maximum fraction of precipitation flowing along the stems
    MOUILLABIL: float = 0.27  # Maximum leaf wetness in mm²/m²
    KSTEMFLOW: float = 0.5  # Not well known, coefficient linking LAI to stemflow

    # Plant height
    HAUTMAX: float = 2.5  # Maximum plant height
    HAUTBASE: float = 0.5  # Base height

    # METHOD_ETP: float = 1  # 1 = crop coefficient / beer law, 2 = energy budget
    CODEBESO: int = 0  # Calculation method: 1 = , 2 = energy budget
    CODEDORMANCE: int = 0  # Dormancy calculation option: 1 = fixed, 2 = Richardson (not implemented), 3 = Bidabe
    CODEINTERCEPT: int = 0  # Interception of water by foliage. 1 = yes, 2 = no

    STADES: dict = field(default_factory=dict)

    S: float = field(init=False)
    FCO2: float = field(init=False)
    FCO2S: float = field(init=False)
    DURVIEI: float = field(init=False)
    DENSITESEM: float = field(init=False)
    ZRAC0: float = field(init=False)


    def __post_init__(self):
        self.S = -np.log(1e-2) / (self.ZLABOUR - self.ZPENTE)
        self.FCO2 = 2 - np.exp(
            np.log(2 - self.ALPHACO2) * (self.CO2 - 350) / (600 - 350)
        )  # Effect of CO2 on RUE, roughly it's 0 if [CO2]=350 and it enhances RUE as [CO2] increases
        self.FCO2S = 1 / (1 + 0.77 * (1 - self.FCO2 / 2.5) * (1 - self.CO2 / 330))
        self.DURVIEI = (
            self.RATIODURVIEI * self.DURVIEF
        )  # Lifespan of early leaves

        sowing_depth = {
            "wheat": 2,
            "vine": 115,
            "corn": 4,
            "tomato": 0.5,
            "rapeseed": 1,
            "pea": 3,
            "strawberry": 1,
            "potato": 12,
            "sorghum": 3,
            "beet": 3,
            "lettuce": 0.5,
        }  # Sowing depth: are these values for sowing depth or root depth?
        root_depth = {
            "wheat": 2,
            "vine": 115,
            "corn": 4,
            "tomato": 0.5,
            "rapeseed": 1,
            "pea": 3,
            "strawberry": 1,
            "potato": 12,
            "sorghum": 3,
            "beet": 3,
            "lettuce": 0.5,
        }  # Root depth on day 0
        # Be careful of the confusion between sowing depth (sowing_depth) and maximum root depth on day 0 (ZRAC0)
        density = {
            "wheat": 300,
            "corn": 7000,
            "vine": 0.4,
            "tomato": 225,
            "rapeseed": 40,
            "pea": 80,
            "strawberry": 15,
            "potato": 10,
            "sorghum": 32000,
            "beet": 10000,
            "lettuce": 30000,
        }  # Seeding density
        perennial_code = {
            "wheat": 1,
            "corn": 1,
            "vine": 2,
            "tomato": 1,
            "rapeseed": 1,
            "pea": 1,
            "strawberry": 1,
            "potato": 1,
            "sorghum": 1,
            "beet": 1,
            "lettuce": 1,
        }
        indeterminate_code = {
            "wheat": 1,
            "corn": 1,
            "vine": 2,
            "tomato": 2,
            "rapeseed": 1,
            "pea": 1,
            "strawberry": 1,
            "potato": 1,
            "sorghum": 1,
            "beet": 2,
            "lettuce": 1,
        }
        self.DENSITESEM = density[self.species]
        self.ZRAC0 = root_depth[self.species]
        self.CODEPERENNE = perennial_code[self.species]
        self.CODEINDETERMIN = indeterminate_code[self.species]

        if self.CODEPERENNE == 1:  # annual plants
            stades_pheno = {
                "0": "germination/levee (ILEV)",
                "1": "dvlp feuilles",
                "2": "tallage",
                "3": "elongation tige princip",
                "4": "gonflement epi",
                "5": "epiaison",
                "6": "floraison",
                "7": "dvlp graines (IDRP)",
                "8": "maturation graines",
                "9": "senescence",
            }
        elif self.CODEPERENNE == 2:  # perennial plants
            stades_pheno = {
                "0": "debourrement (ILEV)",
                "1": "dvlp feuilles",
                "5": "apparition inflo",
                "6": "floraison",
                "7": "dvlp fruits (IDRP)",
                "8": "maturation baies",
                "9": "senescence",
            }


        for key in stades_pheno.keys():
            self.STADES[key] = stades_pheno[key]

    def ini_var(self, dictionary):
        self.__dict__.update(dictionary)


@dataclass
class SoilParams:
    """Soil Parameters."""

    # Soil Profile Description
    EPC_1: int = 30  # Thickness H1
    EPC_2: int = 150  # Thickness H2

    # Soil Water Content
    HMINF_1: float = 0.13  # Humidity at wilting point of Hor.1 [% vol]
    HCCF_1: float = 0.30  # Humidity at field capacity of Hor.1 [cm] or mm of water per cm of soil
    HMINF_2: float = 0.19  # Humidity at wilting point of Hor.2 [cm]
    HCCF_2: float = 0.35  # Humidity at field capacity of Hor.2

    # Soil Properties Driving Evaporation (1st Horizon)
    ARGI: float = (
        8  # Clay content of the 1st layer (%) after decarbonation = used to estimate residual humidity
    )
    DAF_1: float = 1.3  # Apparent density (uncompressed) of the 1st layer (g/cm³)
    Q0: float = 7.5  # [mm] between 0 and 30 = cumulative soil evaporation threshold above which the evaporation rate is reduced = threshold between phase 1 and 2

    # Soil initial state
    HUR1_Init: float = 0.25  # Initial state of water content - Hor.1 [%vol] = proportion of the soil volume
    HUR2_Init: float = 0.30  # Initial state of water content - Hor.2 [%vol]

    ALBEDO: float = 0.2  # NO XML VALUE, albedo of dry and bare soil

    HA: float = field(init=False)
    AEVAP: float = field(init=False)
    ProfSoil: float = field(init=False)

    Z0SOLNU: float = 0.01  # Length of roughness of bare soil. From sols.xml file


    def __post_init__(self):
        self.HA = (
            self.ARGI / 100 / 15 * self.DAF_1
        )  # Residual moisture (mm water per cm of soil) = water content of a soil considered dry = % of clay * apparent density of the surface layer (% vol between [0.001-0.1]). Can we measure it ourselves?
        self.ProfSoil = self.EPC_1 + self.EPC_2  # Maximum soil profile depth [cm]"

    def ini_var(self, dictionary):
        self.__dict__.update(dictionary)


@dataclass
class Constants:
    """Constants"""
    # COEFB: To maintain consistency with ebmax (in g/MJ), we also divide COEFB (in g/MJ) by 100
    COEFB: float = (
        0.0815 / 100
    )  # Value in config > param_gen, saturation effect of radiation on biomass conversion efficiency, in g/MJ
    PROPRAC: float = (
        0.2  # Root mass / above-ground mass at harvest, NO value for Merlot
    )
    Y0MSRAC: float = 0.7  # Minimum root mass at harvest (when above-ground mass is zero)
    LVOPT: float = 0.5  # Root density above which N and H2O absorption are maximum and independent of density = value of the sigmoid plateau
    PARSURRG: float = 0.45  # Ratio of PAR (Photosynthetically Active Radiation) to total radiation

    
    def ini_var(self, dictionary):
        self.__dict__.update(dictionary)


@dataclass
class UserOptions:
    FLAG_HA: bool = (
        False  # False = use the residual water content, True = use wilting point to stop evaporation
    )
    IFINDORM: int = 0  # Date of dormancy break used when CODEDORMANCE = 1. How to implement it to occur at the right moment based on the simulation start date?


@dataclass
class ManageParams:
    ''' Crop management parameters '''
    # Sowing
    IPLT0: int = 295  # sowing day of year. Start of simulation

@dataclass
class StationParams:
    ALBVEG: float = 0.23  # Vegetation albedo, in default > climaisj: always 0.23 in the _sta files, regardless of the crop.
    ACLIM: float = 14  # Wind-related parameter controlling evaporation [mm].
    AANGST: float = 0.18  # Parameter for the Angstrom equation (eq 9.18).
    BANGST: float = 0.62  # Parameter for the Angstrom equation (eq 9.18).
    CODEETP: int = 2  # 1 = ETP provided by external source, 2/3/4 = calculated --> 2 = Penman, 3 = Shuttle, 4 = Pries & Taylor.
    CORECTROSEE: float = 1.0  # Temperature to subtract from Tmin to estimate the dew point.
    LATITUDE: float = 1.0  # Latitude of the station (for calculating photoperiod).
    ZR: float = 129  # Reference height for meteorological data measurements, not provided in the STICS example.

