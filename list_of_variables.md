This document lists all variables used in pySTICS.

# List of parameters
- aangst : float
    > Coefficient of the Angstrom relationship for extraterrestrial radiation. Unit : /
- abscission : float
    > Fraction of senescent leaves falling to the soil.
- aclim : float
    > Climatic component to calculate actual soil evaporation (Brisson & Perrier, 1991). Unit : mm
- adens : float
    > Interplant competition parameter.
- afpf : float
    > Parameter of the logistic function defining sink strength of fruits (indeterminate growth) : relative fruit age at which growth is maximal.
- afruitpot : float
    > Maximal number of set fruits per inflorescence and per degree day (indeterminate growth). Unit : fruit-1.degree days-1
- albedo : float
    > Albedo of the bare dry soil. Unit : /
- albveg : float
    > Albedo of the vegetation. Unit : /
- allocfrmax : float
    > Maximal daily allocation to fruits.
- alphaco2 : float
    > Coefficient accounting for the modification of radiation use efficiency in case of atmospheric CO2 increase.
- alphapt : float
    > Parameter of Priestley-Taylor formula. Unit : /
- ampfroid : float
    > Semi thermal amplitude for vernalising effect. Unit : °C
- argi: float
    > Clay content after decarbonation. Unit : % (0 to 100).
- bangst : float 
    > Coefficient of the Angstrom relationship for extraterrestrial radiation. Unit : /
- bdens : float
    > Minimal plant density above which interplant competition starts. Unit : pl.m-2
- belong : float
    > Parameter of the curve of coleoptile elongation. Unit : degree days-1
- bfpf : float
    > Parameter of the logistic curve defining sink strength of fruits (indeterminate growth): maximum growth rate relative to maximum fruit weight.
- clay : float
    > Soil clay content. Unit : % (0 to 100).
- celong : float
    > Parameter of the plantlet elongation curve.
- cfes : float
    > Parameter defining the soil contribution to evaporation versus depth.
- cfpf : float
    > Parameter of the first potential growth phase of fruit, corresponding to an exponential type function describing the cell division phase.
- cgrain : float
    > Slope of the relationship between grain number and growth rate. Unit : t-1.m2
- cgrainv0 : float
    > Fraction of the maximal number of grains when growth rate is zero.
- codcalinflo : int
    > Option to calculate the inflorescences number: 1 = read in param.par, 2 = calculated at the amf stage.
- codeaumin : int
    > Option to activate the harvest according to grain/fruit water content: 1 = water content > minimum threshold, 2 = water content < maximum threshold.
- codebeso : int
    > Option to calculate water requirements: 1 = k.ETP approach, 2= resistive method.
- codebfroid : int
    > Option to calculate chilling requirements: 1 = no need, 2 = vernalising days (JVC), 3 = development stage (dormancy break and post-dormancy until budding).
- codecalirrig : int
    > Option to activate the automatic calculation of irrigation requirements: 1 = yes, 2 = no.
- codedateapph2o : int
    > Option to calculate irrigation dates according to sum of temperatures: 1 = yes, 2 = no.
- codedate_irrigauto : int
    > Option to activate the beginning and the ending dates in case of automatic irrigation: 1 = dates, 2= crop stages, 3 = nothing
- codedormance : int
    > Option to calculate dormancy and chilling requirements: 1 = forcing, 2 = Richardson, 3 = Bidabe.
- codeetp : int
    > Option to calculate PET: 1 = forced Penman, 2 = calculated Penman, 3= Shuttleworth & Wallace, 4 = Priestley & Taylor.
- codefauche : int
    > Option to activate cuts of forage crops: 1 = yes, 2 = no. Unit : /
- codehypo : int
    > Option to simulate plant emergency: 1 = phase of hypocotyl growth (sown crops), 2 = plantation of plantlets.
- codeintercept : int
    > Option to simulate rainfall interception by leaves: 1 = yes, 2 = no.
- codegermin : int
    > Option to simulate germination: 1 = germination phase, 2 = immediate germination.
- codeindetermin : int
    > Option to simulate the type of leaf growth and fruit growth: 1 = determinate, 2 = undeterminate.
- codeir : int
    > Option to calculate the ratio grain weight/total biomass: 1 = proportional to time, 2 = proportional to thermal time.
- codemodfauche : int
    > Option to define the cutting mode: 1 = automatic calculation depending on phenologic and trophic state, 2 = pre-established calendar in days, 3 = pre-established calendar in degree-days. Unit : /
- codeperenne : int
    > Option to define the crop perenniality: 1 = annual crop, 2 = perennial crop.
- codephot : int
    > Option to define plant photoperiodism: 1 = yes, 2 = no.
- coderetflo : int
    > Option to activate the effect of water stress on development before the stage DRP (filling of harvested organs): 1 = yes, 2 = no.
- codernet : int
    > Option to calculate net radiation: 1 = Brunt method, 2 = Cellier method. 
- codetemp : int
    > Option to calculate thermal time for plant growth: 1 = based on air temperature, 2 = based on crop temperature.
- codetemprac : int
    > Option to calculate thermal time for root growth: 1 = crop temperature, 2 = soil temperature.
- codetremp : int
    > Option to activate heat effect on grain filling: 1 = yes, 2 = no.
- codgelflo : int
    > Option to activate the frost effect at anthesis: 1 = no, 2 = yes.
- codgeljuv : int
    > Option to activate the frost effect on LAI at the juvenile stage: 1 = no, 2 = yes.
- codgellev : int
    > Option to activate the frost effect on plantlet growth: 1 = no, 2 = yes.
- codgelveg : int
    > Option to activate the frost effect on LAI at adult stage: 1 = no, 2 = yes.
- codlocirrig : int
    > Option to define localized irrigation: 1= above the foliage, 2= below the foliage above the soil, 3 = in the soil.
- codrecolte : int
    > Option to define harvest strategy: 1 = at physiological maturity, 2 = according to water content, 3 = according to sugar content, 4 = according to nitrogen content, 5 = according to oil content.
- coefb : float
    > Parameter defining the radiation saturation effect on biomass conversion efficiency. Unit : g.MJ-1
- coefmshaut : float
    > Ratio of crop biomass to useful cutting height of crops. Unit : t.ha-1.m-1
- contrdamax : float
    > Maximal reduction factor applied to root growth rate due to soil strengthness (high bulk density).
- corectrosee : float
    > Temperature to substract to Tmin to estimate dew point temperature (in case of missing air humidity data). Unit : °C
- co2 : float
    > Atmospheric CO2 concentration. Unit : ppm
- croirac : float
    > Elongation rate of the root apex. Unit : cm.degree days-1
- daf : numpy array
    > Bulk density of fine earth fraction for each soil cm. Unit : g.cm-3
- daf_1 : float
    > Bulk density of fine earth fraction in soil layer 1. Unit : g.cm-3
- daf_2 : float
    > Bulk density of fine earth fraction in soil layer 2. Unit : g.cm-3
- daf_3 : float
    > Bulk density of fine earth fraction in soil layer 3. Unit : g.cm-3
- daf_4 : float
    > Bulk density of fine earth fraction in soil layer 4. Unit : g.cm-3
- daf_5 : float
    > Bulk density of fine earth fraction in soil layer 5. Unit : g.cm-3
- daseuilbas : float
    > Bulk density of soil above which root growth is maximal. Unit : g.cm-3
- daseuilhaut : float
    > Bulk density of soil above which root growth becomes impossible. Unit : g.cm-3
- datedeb_irrigauto : int
    > Starting date of automatic irrigations. Unit : julian day
- datefin_irrigauto : int
    > Ending date of automatic irrigations. Unit : julian day
- densinitial_1 : float
    > Initial root density of soil layer 1. Unit : cm.cm-3
- densinitial_2 : float
    > Initial root density of soil layer 2. Unit : cm.cm-3
- densinitial_3 : float
    > Initial root density of soil layer 3. Unit : cm.cm-3
- densinitial_4 : float
    > Initial root density of soil layer 4. Unit : cm.cm-3
- densinitial_5 : float
    > Initial root density of soil layer 5. Unit : cm.cm-3
- densitesem : int
    > Plant sowing density. Unit : pl.m-2
- depth : int
    > Soil depth. Unit : cm
- deshydbase : float
    > Rate of change of water content in fruits (FW) vs thermal time (>0 or <0). Unit : g.g-1.degree days-1
- diftherm : float
    > Soil thermal diffusivity. Unit : cm2.s-1
- dlaimaxbrut : float
    > Maximum rate of net daily increase of LAI. Unit : m2.pl-1.degree days-1
- dlaimin : float
    > Accelerating parameter for the lai growth rate.
- doseirrigmin : int
    > Minimal amount of daily irrigation. Unit : mm
- dosimx : int
    > Maximum amount of irrigation water applied daily (mode automatic irrigation). Unit : mm
- dureefruit : int
    > Duration of the fruit between onset and physiological maturity. Unit : degree days
- durvief : int
    > Maximal lifespan of an adult leaf. Unit : Q10
- efcroijuv : float
    > Maximum radiation use efficiency during the juvenile phase (LEV-AMF). Unit : g.MJ-1
- efcroirepro : float
    > Maximum radiation use efficiency during the grain filling phase (DRP-MAT). Unit : g.MJ-1
- efcroiveg : float
    > Maximum radiation use efficiency during the vegetative stage (AMF-DRP). Unit : g.MJ-1
- effirr : float
    > Irrigation efficiency.
- efremobil : float
    > Efficiency of use of carbohydrates in storage organs of perennials.
- elmax : float
    > Maximum elongation of the coleoptile in darkness condition. Unit : cm
- epc_1 : float
    > Thickness of soil layer 1. Unit : cm
- epc_2 : float
    > Thickness of soil layer 2. Unit : cm
- epc_3 : float
    > Thickness of soil layer 3. Unit : cm
- epc_4 : float
    > Thickness of soil layer 4. Unit : cm
- epc_5 : float
    > Thickness of soil layer 5. Unit : cm
- extin : float
    > Extinction coefficient of photosynthetic active radiation in the canopy. Unit : /
- forage : boolean
    > Code which determines if the plant is a forage crop or not. Unit : /
- gamma : float
    > Psychrometrc constant. Unit : mb.°C-1
- hautbase : float
    > Basal height of crop. Unit : m
- hautcoupe : float
    > Cut height for forage crops (calendar fixed). Unit : m
- hautmax : float
    > Maximum height of crop. Unit : m
- hccf_1 : float
    > Gravimetric water content at field capacity for soil layer 1. Unit : % TODO : checker %vol ou %pond
- hccf_2 : float
    > Gravimetric water content at field capacity for soil layer 2. Unit : % TODO : checker %vol ou %pond
- hccf_3 : float
    > Gravimetric water content at field capacity for soil layer 3. Unit : % TODO : checker %vol ou %pond
- hccf_4 : float
    > Gravimetric water content at field capacity for soil layer 4. Unit : % TODO : checker %vol ou %pond
- hccf_5 : float
    > Gravimetric water content at field capacity for soil layer 5. Unit : % TODO : checker %vol ou %pond
- herbaceous : bool
    > Code to define if the crop is herbaceous or ligneous. True if herbaceous, False if ligneous.
- hinitf_1 : float
    > Initial gravimetric water content of soil layer 1. Unit : mm water.cm soil-1
- hinitf_2 : float
    > Initial gravimetric water content of soil layer 2. Unit : mm water.cm soil-1
- hinitf_3 : float
    > Initial gravimetric water content of soil layer 3. Unit : mm water.cm soil-1
- hinitf_4 : float
    > Initial gravimetric water content of soil layer 4. Unit : mm water.cm soil-1
- hinitf_5 : float
    > Initial gravimetric water content of soil layer 5. Unit : mm water.cm soil-1
- hminf_1 : float
    > Gravimetric water content at wilting point for soil layer 1. Unit : % TODO : checker %vol ou %pond
- hminf_2 : float
    > Gravimetric water content at wilting point for soil layer 2. Unit : % TODO : checker %vol ou %pond
- hminf_3 : float
    > Gravimetric water content at wilting point for soil layer 3. Unit : % TODO : checker %vol ou %pond
- hminf_4 : float
    > Gravimetric water content at wilting point for soil layer 4. Unit : % TODO : checker %vol ou %pond
- hminf_5 : float
    > Gravimetric water content at wilting point for soil layer 5. Unit : % TODO : checker %vol ou %pond
- h2ofrvert : float
    > Water content of fruits before the beginning of dehydration (FW). Unit : g.g-1
- h2ograinmax : float
    > Maximal water content of fruits at harvest (FW). Unit : g.g-1
- h2ograinmin : float
    > Minimal water content of fruits at harvest (FW). Unit : g.g-1
- iplt0 : int
    > Date of sowing. Unit : julian day
- irmax : float
    > Maximum harvest index.
- irrigation_interventions : dict
    > Days and amounts of irrigation when codecalirrig = 2. Units : julian day / mm
- irrlev : float
    > Amount of irrigation applied automatically on the sowing day to allow germination when the model calculates irrigation. Unit : mm
- julfauche : list
    > Julian days for cutting intervention (if codemodfauche = 2). Unit : julian day
- jvc : int
    > Number of vernalising days or dormancy units. Unit : days
- jvc_mini : int
    > Minimum number of vernalising days. Unit : days
- khaut : float
    > Extinction coefficient connecting LAI to crop height. Unit : /
- kmax : float
    > Maximum crop coefficient for water requirements (= MET/PET).
- kstemflow : float
    > Extinction coefficient connecting LAI to stemflow.
- laicomp : float
    > LAI above which competition between plants starts. Unit : m2.m-2
- lai0 : float
    > Initial leaf area index. Unit : m2.m-2
- latitude : float
    > Latitude of the plot. Unit : °
- longitude : float
    > Longitude of the plot. Unit : °
- lvopt : float
    > Root length density (RLD) above which water and N uptake are maximum and independent of RLD. Unit : cm.cm-3
- masec0 : float
    > Initial plant biomass (if the option to simulate N and C reserves is not activated). Unit : t.ha-1
- mouillabil : float
    > Maximum wettability of leaves. Unit : mm.m-1
- nbfgellev : int
    > Leaf number at the end of the juvenile phase (frost sensitivity). Unit : pl-1
- nbgrmax : int
    > Maximum number of fruits per surface area. Unit : m-2
- nbgrmin : int
    > Minimum number of fruits per surface area. Unit : m-2
- nbinflo : int
    > Imposed number of inflorescences per plant. Unit : pl-1
- nbjgerlim : int
    > Maximum number of days after grain imbibition allowing full germination. Unit : days
- nbjgrain : int
    > Number of days used to compute the number of viable grains. Unit : days
- nboite : int
    > Number of boxes or age classes of fruits used to calculate fruit growth for undeterminate crops.
- nlevlim1 : int
    > Number of days after germination after which plant emergence is reduced. Unit : days
- nlevlim2 : int
    >  	Number of days after germination after which plant emergence is impossible. Unit : days
- parsurrg : float
    > Fraction of photosynthetically active radiation in global radiation (PAR/RG). Unit : /
- pentinflores : float
    > Parameter used to calculate the inflorescences number. Unit : kg-1
- pentlaimax : float
    > Parameter of the logistic curve of LAI growth.
- pgrainmaxi : float
    > Maximum grain weight (at 0% water content). Unit : g
- phobase : int
    > Basal photoperiod. Unit : hours
- phosat : int
    > Saturating photoperiod. Unit : hours
- phyllotherme : int
    > Thermal duration between the apparition of two successive leaves on the main stem. Unit : degree days
- potgermi : float
    > Soil water potential below which seed imbibition is impeded. Unit : MPa
- profsem : float
    > Depth of sowing. Unit : cm
- profsoil : int
    > Soil depth. Unit : cm
- propjgermin : float
    > Minimal fraction of the duration nbjgerlim when the temperature is higher than the temperature threshold Tdmax. Unit : 0 to 1.
- proprac : float
    > Ratio of root mass to aerial mass at harvest.
- propres : float
    > Maximal fraction of the biomass reserves that can be mobilized from aerial organs in all crops.
- propresp : float
    > Maximal fraction of the biomass reserves that can be mobilized from storage organs in perennials.
- psihucc : float
    > Soil water potential corresponding to field capacity. Unit : MPa
- psihumin : float
    > Soil water potential corresponding to wilting point. Unit : MPa
- psisto : float
    > Potential of stomatal closing (absolute value). Unit : bars
- psiturg : float
    > Potential of the beginning of decrease of the cellular extension (absolute value). Unit : bars
- q0 : float
    > Cumulative soil evaporation above which evaporation rate is decreased. Unit : mm
- q10 : float
    > Q10 used for the dormancy break calculation.
- rapsenturg : float
    > Threshold soil water content active to simulate water senescence stress as a proportion of the turgor stress.
- ratiodurviei : float
    > Life span of early leaves expressed as a fraction of the life span of the last leaves emitted DURVIEF.
- ratiol : float
    > Water stress index below which irrigation is started in automatic mode (0 in manual mode).
- ratiosen : float
    > Fraction of senescent biomass (relative to total biomass).
- rayon : float
    > Average radius of the roots. Unit : cm
- remobres : float
    > Fraction of daily remobilisable C reserves. Unit : day-1
- resplmax : float
    > Maximal reserve of biomass. Unit : t.ha-1
- rsmin : float
    > Minimal stomatal resistance of leaves. Unit : s.m-1
- sand : float
    > Soil sand content. Unit : % (0 to 100).
- sensiphot : int
    > Index of photoperiod sensitivity (1 = insensitive, 0 = highly sensitive).
- sensrsec : float
    > Index of root sensitivity to drought (1 = insensitive, 0 = highly sensitive).
- silt : float
    > Soil silt content. Unit : % (0 to 100).
- slamax : int
    > Maximum SLA (specific leaf area) of green leaves. Unit : cm2.g-1
- slamin : int
    > Minimum SLA (specific leaf area) of green leaves. Unit : cm2.g-1
- soc : float
    > Amount of soil organic C (= Chumt + Cb) over the depth with an active biological activity. Unit : % (0 to 100).
- species : string
    > Species of the plant simulated.
- spfrmax : float
    > Maximal sources/sinks value allowing the trophic stress calculation for fruit onset.
- spfrmin : float
    > Minimal sources/sinks value allowing the trophic stress calculation for fruit onset.
- stade0 : string
    > Crop stage at the beginning of simulation.
- stamflax : int
    > Cumulative thermal time between the stages AMF (maximum acceleration of leaf growth, end of juvenile phase, 1st npde) and LAX (maximum leaf area index, end of leaf growth, flag-leaf). Unit : degree days
- stdordebour : int
    > Cumulative thermal time between the dormancy break and the bud break. Unit : degree days
- stades : dict
    > List of phenological stages and name.
- stdrpdes : int
    > Cumulative thermal time between the DRP stage (starting date of filling of harvested organs) and DEBDES (date of onset of water dynamics in harvested organs). Unit : degree days
- stdrpmat : int
    > Cumulative thermal time between the stages DRP (starting date of filling of harvested organs) and MAT (maturity). Unit : degree days
- stemflowmax : float
    > Maximal fraction of rainfall flowing down along the stems.
- stflodrp : int
    > Cumulative thermal time between FLO (anthesis) and DRP (starting date of filling of harvested organs) (only for indication). Unit : degree days
- stlaxsen : int
    > Cumulative thermal time between the stages LAX (maximum leaf area index, end of leaf growth, flag-leaf) and SEN (beginning of leaf senescence). Unit : degree days
- stlevamf : int
    > Cumulative thermal time between the stages LEV (emergence) and AMF (maximum acceleration of leaf growth, end of juvenile phase, 1st npde). Unit : degree days
- stlevdrp : int
    > Cumulative thermal time between the stages LEV (emergence) and DRP (starting date of filling of harvested organs). Unit : degree days
- stoprac : string
    > Stage when root growth stops (LAX= maximum leaf area index, end of leaf growth or SEN=beginning of leaf senescence).
- stpltger : int
    > Cumulative thermal time allowing germination. Unit : degree days
- stressdev : float
    > Maximum phasic delay allowed due to stresses.
- stdordebour : int
    > Cumulative thermal time between the dormancy break and the bud break. Unit : degree days
- stdordebour : int
    > Cumulative thermal time between the dormancy break and the bud break. Unit : degree days
- stressdev : float
    > Maximum phasic delay allowed due to stresses.
- swfacmin : float
    > Minimal value for drought stress index (turfac, swfac, senfac).
- tcmax : float
    > Maximum temperature at which growth ceases. Unit : °C
- tcmin : float
    > Minimum temperature at which growth ceases. Unit : °C
- tcxstop : float
    > Temperature beyond which foliar growth stops. Unit : °C
- tdebgel : float
    > Temperature below which frost affects plant growth. Unit : °C
- tdmax : float
    > Maximum temperature above which development stops. Unit : °C
- tdmin : float
    > Minimum temperature below which development stops. Unit : °C
- tdmaxdeb : float
    > Maximal thermal threshold for hourly calculation of phasic duration between dormancy and bud breaks. Unit : °C
- tdmindeb : float
    > Minimal thermal threshold for hourly calculation of phasic duration between dormancy and bud breaks. Unit : °C
- temax : float
    > Maximal temperature above which plant growth stops. Unit : °C
- temin : float
    > Minimum temperature for development. Unit : °C
- tempdeshyd : float
    > Increase in fruit dehydration rate due to the increase in crop temperature (Tcult-Tair). Unit : %.°C-1
- tempfauche : list
    > Cumulative thermal time between two cuts of forage crops. Unit : degree days
- teopt : float
    > Optimal temperature (1/2) for plant growth. Unit : °C
- teoptbis : float
    > Optimal temperature (2/2) for plant growth. Unit : °C
- tfroid : float
    > Optimal temperature for vernalisation. Unit : °C
- tgelflo10 : float
    > Temperature resulting in 10% of frost damages on flowers and fruits. Unit : °C
- tgelflo90 : float
    > Temperature resulting in 90% of frost damages on flowers and fruits. Unit : °C
- tgeljuv10 : float
    > Temperature resulting in 10% of frost damage on LAI (juvenile stage). Unit : °C
- tgeljuv90 : float
    > Temperature resulting in 90% of frost damage on LAI (juvenile stage). Unit : °C
- tgellev10 : float
    > Temperature resulting in 10% of frost damages on plantlet. Unit : °C
- tgellev90 : float
    > Temperature resulting in 90% of frost damages on plantlet. Unit : °C
- tgelveg10 : float
    > Temperature resulting in 10% of frost damage on LAI (adult stage). Unit : °C
- tgelveg90 : float
    > Temperature resulting in 90% of frost damage on LAI (adult stage). Unit : °C
- tgmin : float
    > Minimum temperature below which emergence is stopped. Unit : °C
- tigefeuil : float
    > Ratio stem (structural part)/leaf.
- tletale : float
    > Lethal temperature for the plant. Unit : °C
- tmaxrempl : float
    > Maximal temperature above which grain filling stops. Unit : °C
- tminremp : float
    > Minimal temperature below which grain filling stops. Unit : °C
- udlaimax : float
    > Ulai from which the rate of leaf growth decreases.
- vitircarb : float
    > Rate of increase of the C harvest index vs time. Unit : g.g-1.day-1
- vitircarbt : float
    > Rate of increase of the C harvest index vs thermal time. Unit : g.g-1.degree day-1
- vlaimax : float
    > Ulai at the inflexion point of the function DELTAI=f(ULAI).
- variety : string
    > Variety of the plant simulated.
- y0msrac : float
    > Minimal amount of root mass at harvest (when aerial biomass is nil). Unit : t.ha-1
- zesx : int
    > Maximal soil depth affected by soil evaporation. Unit : cm
- zlabour : int
    > Depth of ploughing (reference profile). Unit : cm
- zpente : int
    > Depth at which root density is 50% of the surface root density (reference profile). Unit : cm
- zprlim : int
    > Maximum depth of the root profile (reference profile). Unit : cm
- zr : int
    > Reference height of meteorological data measurement. Unit : m
- zrac0 : float
    > Initial depth of root apex of the crop. Unit : cm
- z0solnu : float
    > Roughness length of bare soil. Unit : m   


# List of weather inputs variables
- temp_max : 
    > Daily maximum air temperature. Unit : °C
- temp_min : 
    > Daily minimum air temperature. Unit : °C
- trg : float
    > Active radiation (entered or calculated). Unit : MJ.m-2
    
# List of computed variables
- aevap : float
   	> Soil evaporation parameter combining climatic and soils aspects. Unit : mm
- airg : float
    > Daily amount of irrigation water. Unit : mm
- albedolai : float
    > Albedo of the crop including soil and vegetation. Unit : /
- albsol : float
    > Albedo of the soil. Unit : /
- albsolhum : float
    > Modified soil_albedo to take into account humidity. Unit : /
- amf : list
    > End of juvenile phase indexes : 1 = end of juvenile phase has been reached, 0 = end of juvenile phase has not been reached. Unit : /
- amplsurf : float
    > Daily thermal amplitude at the surface. Unit : °C
- coeflev : float
    > Ratio between the emerged and the germinated density depending on non-optimal water content and temperature conditions. Unit : /
- cu : float
    > Degree days based on Q10 and air temperature. Unit : °C.d
- cumlracz : float
    > Cumulative length of active roots per soil surface. Unit : cm.cm-2
- dayLAIcreation : int
    > Creation day of oldest non senescent deltai. Unit : julian day
- debdes : int
    > Onset of water dynamics indexes : 1 = onset of water dynamics has been reached, 0 = onset of water dynamics has not been reached. Unit : /
- deltai : float
   	> Daily increase of the leaf area index. Unit : m2.m-2
- deltai_dens : float
    > Density component of deltai. Unit : pl.m-2
- deltai_dev : float
   	> Phasic development component of deltai. Unit : m2.degree day-1
- deltai_stress : float
    > Stress component of deltai. Unit : /
- deltai_t : float
    > Thermal component of deltai. Unit : degree days
- deltaz : float
    > Root front growth. Unit : cm
- deltaz_stress : float
   	> Stress component of deltaz. Unit : /
- deltaz_t : float
   	> Thermal component of deltaz. Unit : degree days
- densite : float
    > Actual plant density. Unit : pl.m-2
- dh : float
    > Displacement height. Unit : m
- dltags : float
   	> Daily growth rate of the plantlets. Unit: t.ha-1
- dltaisen : float
    > Daily change in the senescent leaf area index. Unit : m2.m-2
- dltams : float
    > Daily growth rate of the plant. Unit : t.ha-1
- dltamsen : float
    > Daily senescence rate of the plant. Unit : t.ha-1
- dltamstombe : float
    > Daily fallen biomass. Unit : t.ha-1
- drain : float
    > Daily amount of water drained at the base of the soil profile. Unit : mm
- drp : list
    > Fruit filling indexes : 1 = fruit filling onset has been reached, 0 = fruit filling onset has not been reached. Unit : /
- durage : float
    > Natural lifespan of leaves. Unit : °C
- durvie : float
   	> Actual life span of the leaf surface. Unit : °C
- durviei : float
    > Lifespan of a young leal (at the AMF stage) expressed as a proportion of DURVIEF.
- ebmax : float
   	> Maximum value of radiation use efficiency. Unit : cg.MJ-1
- edirect : float
   	> Actual evaporation of soil water and water intercepted by leaves. Unit : mm
- edirectm : float
   	> Maximal evaporation of soil water and water intercepted by leaves. Unit : mm
- efda : float
    > Reduction factor on root growth due to physical constraint (through bulk density). Unit : /
- efdensite : float
   	> Density factor on leaf area growth. Unit : /
- elong : float
    > Increase in fruit dehydration rate due to the increase in crop temperature (Tcult-Tair). Unit : cm
- emd : float
    > Daily amount of water directly evaporated after leaves interception. Unit : mm
emissa : float
- eo : float
   	> Intermediary variable for the computation of evapotranspiration. Unit : mm
- eop : float
   	> Daily maximum transpiration flux. Unit : mm
- eos : float
   	> Daily maximum evaporation flux. Unit : mm
- ep : float
   	> Daily actual transpiration flux. Unit : mm
- esol : float
   	> Daily actual soil evaporation flux. Unit : mm
- et : float
    > Actual evapotranspiration. Unit : mm
- etp : float
    > Daily potential evapotranspiration. Unit : mm
- fco2 : float
    > Specie-dependant CO2 effect on radiation use efficiency.
- fco2s : float
    > Specie-dependant CO2 effect onstomate closure.
- fgelflo : float
    > Reduction factor on the number of fruits due to frost. Unit : /
- fgellev : float
    > Frost index acting on plant density during the plantlet phase. Unit : /
- findorm : int
    > Dormancy break indexes : 1 = dormancy break has been reached, 0 = dormancy break has not been reached. Unit : /
- flagrain : int
    > Index of raining day (1 = raining day, 0 = no rain). Unit : /
- flo : list
    > Flowering indexes : 1 = flowering has been reached, 0 = flowering has not been reached. Unit : /
- fracinsol : float
    > Insolation fraction. Unit : 0 to 1
- fstressgel : float
   	> Reduction factor on leaf growth due to frost. Unit : /
- ftemp : float
    > Thermal reduction factor on biomass growth. Unit : /
- ftempremp : float
    > Thermal stress index as a function of temperature using cardinal temperatures (tminremp et tmaxremp). Unit : 0 to 1
- gdh : float
    > Growing degree hours. Unit : °C
- ger : int
    > Germination indexes : 1 = germination has been reached, 0 = germination has not been reached. Unit : /
- ha : float
   	> Residual soil water content in the sed bed. Unit : mm.cm-1
- hauteur : float- 
    > Crop height. Unit : m
- ircarb : float
   	> Carbon harvest index. Unit : /
- jvi : int
    > Vernalizing contribution of a given day. Unit : days
- lai : float
    > Leaf area index. Unit : m2.m-2
- laisen : float
   	> Leaf area index of senescent leaves . Unit : m2.m-2
- lai_prev : float
    > Leaf area index. Unit : /
- lax : list
    > End of leaf growth indexes : 1 = end of leaf growth has been reached, 0 = end of leaf growth has not been reached. Unit : /
- let : int
    > Plantlet stage indexes : 1 = plantlet stage has been reached, 0 = plantlet stage has not been reached. Unit : /
- lev : list
    > Emergence indexes : 1 = emergence has been reached, 0 = emergence has not been reached. Unit : /
- lev_i : int
    > Emergence index of current day : 1 = emergence has been reached, 0 = emergence has not been reached. Unit : /
- lev_i_prev : int
    > Emergence index of previous day : 1 = emergence has been reached on previous day, 0 = emergence has not been reached. Unit : /
- mafeuiljaune : float
    > Biomass of yellow leaves. Unit : t.ha-1 
- mafeuiltombe : float
    > Cumulated fallen biomass. Unit : t.ha-1
- mafeuilverte : float
    > Biomass of green leaves. Unit : t.ha-1
- mafruit : float
    > Biomass of harvested organs. Unit : t.ha-1
- masec : float
   	> Biomass of aboveground plant. Unit : t.ha-1
- masecneo : float
   	> Biomass of newly-formed organs. Unit : t.ha-1
- masectot : float
    > Total plant biomass (aerials + roots + perennial organs). Unit : t.ha-1
- mat : int
    > Maturation indexes : 1 = maturation has been reached, 0 = maturation has not been reached. Unit : /
- moist : int
    > Moistening indexes : 1 = moistening has been reached, 0 = moistening has not been reached. Unit : /
- mouill : float
    > Water retained on the foliage. Unit : mm
- msneojaune : float
   	> Newly-formed senescent biomass. Unit : t.ha-1
- msrec_fou : float
   	> Biomass of harvested forage. Unit : t.ha-1
- msres : float
    > Non-senescent residual dry matter. Unit : t.ha-1
- msresjaune : float
   	> Senescent residual dry matter. Unit : t.ha-1
- nbfeuille : int
    > Number of leaves on main stem. Unit : /
- nbgraingel : float
    > Number of frozen grains per square meter. Unit : grains.m-2
- nbgrains : float
    > Number of grains per square meter. Unit : grains.m-2
- nbjgrauto : float
    > Days of autotrophy for a moistened seed. Unit : days
- nbjhumec : float
    > Maximal period that seedcan be in a moist status without seed death occurs. Unit : days
- nou : list
    > End of setting indexes : 1 = end of setting has been reached, 0 = end of setting has not been reached. Unit : /
- pgrain : float
    > Weight of grains per square meter. Unit : grains.m-2
- pgraingel : float
    > Weight of frozen grains per square meter. Unit : g.m-2
- phoi : float
    > Photoperiod. Unit : hour
- psibase : float
   	> Predawn leaf water potential. Unit : MPa
- psisols : float
    > Retention curve parameter 1. Unit : Mpa
- raint : float
    > Photosynthetic active radiation intercepted by the canopy. Unit : MJ.m-2
- ratm : float
    > Atmospheric radiation. Unit : MJ.m-2
- rfpi : float
   	> Reduction factor on plant development due to photoperiod. Unit: /
- rfvi : float
    > Reduction factor on plant development due to vernalization. Unit : /
- rgex : float
    > Extraterrestrial radiation. Unit : MJ.m-2
- rglo : float
    > Long wave radiation. Unit : MJ.m-2
- rnet : float
    > Net radiation. Unit : MJ.m-2
- senfac : float
   	> Reduction factor on leaf life span due to water stress (increasing senescence rate). Unit : /
- senstress : float
    > Stress index affecting senescence. Unit : 0 to 1
- seos : float
    > Cumulated daily maximum evaporation flux with phase 1/2 reduction. Unit : mm
- somfeuille : float
    > Cumulative thermal time during frost sensitivity period. Unit : °C.d
- somger : float
    > Growing degree-days from planting in the seed bed. Unit : °C
- somsen : float
    > Current thermal time for senescence. Unit : °C
- somtemp : float
    > Sum of temperatures expressed in Q10 =sum (2.0 ** (udevair ou udevcult / 10.)). Unit : °C.d
- stemflow : float
    > Daily amount of water runoff along the stem. Unit : mm
- stopfeuille : string
    > Phenological stage when leaf growth is stopped.
- sum_upvt_post_lev : float
    > Cumulated development unit (0 before emergence). Unit : degree day
- swfac : float
   	> Stomatic water stress index. Unit : /
- tcult : float
    > Crop surface temperature. Unit : °C
- tcultmax : float
    > Crop max surface temperature. Unit : °C
- tcult_tmp : float
    > Temporary crop surface temperature during iterative calculation. Unit : °C
- teaugrain : float
    > Water content of (harvested) organs. Unit : g.g-1
- temp : float
    > Daily average air temperature ((min+max)/2). Unit : °C
- teta : float
    > Available water content in the root zone. Unit : cm3.cm-3
- tetsen : float
    > Threshold soil water content accelering senescence. Unit : mm.cm-1
- tetstomate : float
   	> Threshold of soil water content limiting transpiration and photosynthesis. Unit : mm.cm-1
- teturg : float
   	> Threshold of soil water content limiting the growth of leaves (in surface area). Unit : mm.cm-1
- tpm : float
    > Water vapour pressure in air. Unit : hPa
- trg : float
    > Active radiation (entered or calculated). Unit : MJ.m-2
- trr : float
   	> Daily rainfall. Unit : mm
- turfac : float
    > Turgescence water stress index. Unit : /
- tvent : float
   	> Mean daily wind speed at 2 m high above soil. Unit : /
- udevcult : float
   	> Effective temperature for crop development, computed with tcult. Unit : degree day
- udevcult_lev : float
   	> Effective temperature during frost sensitivity period, computed with tcult. Unit : degree day
- ulai : float
    > Relative development unit for leaf growth. Unit : /
- upvt_post_lev : float
   	> Development unit (0 before emergence). Unit : degree day
- vitmoy : float
   	> Mean canopy growth rate. Unit : g.m-2
- zdemi : float
    > Root depth that ensures at least an extraction near the soil surface of 20% of the water available. Unit : cm
- znonli : float
    > Root depth if no obstacle nor phenological stop. Unit : cm
- zrac : float
    > Root front depth. Unit : cm
- z0 : float
    > Crop roughness. Unit : m


# List of data class objects
- crop : CropParams object
    > Object containing the parameters for the crop being simulated.
- soil : SoilParams object
    > Object containing the parameters for the soil being used in the simulation.
- constants : Constants object
    > Object containing the constants parameters of the simulation.
- initial : Initial object
    > Object containing the initial parameters for the simulation.
- station : Station object
    > Object containing the station parameters for the simulation.


# List of output matrixes
- outputs : pandas DataFrame
    > Output dataframe containing all computed variables for each day.
- lracz : numpy array
    > Array containing root length density for each soil depth (columns) for each day (lines)
- hur : numpy array
    > Array containing soil water content for each soil depth (columns) for each day (lines)
- wi : numpy array
    > Array containing soil water availability index for each soil depth (columns) for each day (lines)
- esz : numpy array
    > Array containing soil evaporation contribution for each soil depth (columns) for each day (lines)
- tsol : numpy array
    > Array containing soil temperature for each soil depth (columns) for each day (lines)
- humirac : numpy array
    > Array containing water deficit index (affecting emergence and root growth) for each soil depth (columns) for each day (lines)
- humpotsol : numpy array
    > Array containing soil water content (associated to water potential) for each soil depth (columns) for each day (lines)
- psisol : numpy array
    > Array containing soil potential for each soil depth (columns) for each day (lines)
- racinepsi : numpy array
    > Array containing root length density participating in predawn potential (in moist layers) for each soil depth (columns) for each day (lines)
- amplz : numpy array
    > Array containing thermal amplitude for each soil depth (columns) for each day (lines)