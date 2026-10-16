# has NH3 in it to see this interaction
# Ni is the only surface with ammonia library

# Data sources
database(
    thermoLibraries=['surfaceThermoNi111', 'primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC','DFT_QCI_thermo', 'auto'], 
    reactionLibraries = [ ('Surface/Methane/Deutschmann_Ni', True), ('Surface/Ammonia/Duan_Ni111', True),('O2_dissociative_adsorption', False), "auto"], 
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['surface','default', "auto"],
    kineticsEstimator = 'rate rules',
)

catalystProperties(metal="Ni111")

# List of species

forbidden(
    label = "surface_singlet_carbene",
    structure = adjacencyListGroup(
        """
        1 C u0 p1 c0
        2 X u0 p0 c0
        """
    ),
)

species(
    label='H2',
    reactive=True,
    structure=SMILES("[HH]"),
)

species(
    label='H2O',
    reactive=True,
    structure=SMILES("O"),
)

species(
    label='CO2',
    reactive=True,
    structure=SMILES("C(=O)=O"),
)

species(
    label='CO',
    reactive=True,
    structure=SMILES("[C-]#[O+]"),
)

species(
    label='CH4',
    reactive=True,
    structure=SMILES("[CH4]"),
)

species(
   label='O2',
   reactive=True,
   structure=SMILES("[O][O]"),
)

species(
    label='N2',
    reactive=True,
    structure=SMILES("N#N"),
)

species(
    label='vacantX',
    reactive=True,
    structure=adjacencyList("1 X u0"),
)

species(
    label='NH3',
    reactive=True,
    structure=SMILES("N"),
)
species(
    label='N2O',
    reactive=True,
    structure=SMILES("[O-][N+]#N"),
)

species(
    label='C2H6',
    reactive=True,
    structure=SMILES("CC"),
)

species(
    label='C2H4',
    reactive=True,
    structure=SMILES("C=C"),
)


species(
    label='C2H2',
    reactive=True,
    structure=SMILES("C#C"),
)

### Pressure 1

# Reaction systems
surfaceReactor(
    # range of temps
    temperature=[(650,'K'),(1200,'K')],
    # only works with one pressure, more stuff happening when higher
    initialPressure=(1.0, 'bar'),
    # number of simulations
    nSims = 5,
    # gas mole fractions, will be normalized
    # adding CO2 and H2 on their own 
    initialGasMoleFractions = {"H2O": 10,
                               "CO2": 0.05,
                               "H2": 190,
                               "CO": 1,
                               "CH4": 0.8,
                               "N2": 0.85,
                               "O2": 0.15,
                               "NH3": 0.1,},
    # surface coverage, this is proportion, so this is wide open
    initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "H2":0.99,},
    terminationTime=(1000, 's'),
    #terminationRateRatio=0.01,
)
surfaceReactor(
    # range of temps
    temperature=[(650,'K'),(1200,'K')],
    # only works with one pressure, more stuff happening when higher
    initialPressure=(1.0, 'bar'),
    # number of simulations
    nSims = 5,
    # gas mole fractions, will be normalized
    # adding CO2 and H2 on their own 
    initialGasMoleFractions = {"H2O": 0,
                               "CO2": 0,
                               "H2": 8,
                               "CO": 2,
                               "CH4": 0,
                               "N2": 0,
                               "O2": 0,},
    # surface coverage, this is proportion, so this is wide open
    initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "H2":0.99,},
    terminationTime=(1000, 's'),
    #terminationRateRatio=0.01,
)
surfaceReactor(
    # range of temps
    temperature=[(650,'K'),(1200,'K')],
    # only works with one pressure, more stuff happening when higher
    initialPressure=(1.0, 'bar'),
    # number of simulations
    nSims = 5,
    # gas mole fractions, will be normalized
    # adding CO2 and H2 on their own 
    initialGasMoleFractions = {"H2O": 0,
                               "CO2": 2,
                               "H2": 8,
                               "CO": 0,
                               "CH4": 0,
                               "N2": 0,
                               "O2": 0,},
    # surface coverage, this is proportion, so this is wide open
    initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "H2":0.99,},
    terminationTime=(1000, 's'),
    #terminationRateRatio=0.01,
)
surfaceReactor(
    # range of temps
    temperature=[(650,'K'),(1200,'K')],
    # only works with one pressure, more stuff happening when higher
    initialPressure=(1.0, 'bar'),
    # number of simulations
    nSims = 5,
    # gas mole fractions, will be normalized
    # adding CO2 and H2 on their own 
    initialGasMoleFractions = {"H2O": 0,
                               "CO2": 0,
                               "H2": 9,
                               "CO": 0,
                               "CH4": 1,
                               "N2": 0,
                               "O2": 0,},
    # surface coverage, this is proportion, so this is wide open
    initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "H2":0.99,},
    terminationTime=(1000, 's'),
    #terminationRateRatio=0.01,
)
surfaceReactor(
    # range of temps
    temperature=[(650,'K'),(1200,'K')],
    # only works with one pressure, more stuff happening when higher
    initialPressure=(1.0, 'bar'),
    # number of simulations
    nSims = 5,
    # gas mole fractions, will be normalized
    # adding CO2 and H2 on their own 
    initialGasMoleFractions = {"H2O": 0,
                               "CO2": 0,
                               "H2": 9,
                               "CO": 0,
                               "NH3": 1,
                               "N2": 0,
                               "O2": 0,},
    # surface coverage, this is proportion, so this is wide open
    initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "H2":0.99,},
    terminationTime=(1000, 's'),
    #terminationRateRatio=0.01,
)


### Pressure 2

# Reaction systems
surfaceReactor(
    # range of temps
    temperature=[(650,'K'),(1200,'K')],
    # only works with one pressure, more stuff happening when higher
    initialPressure=(100.0, 'bar'),
    # number of simulations
    nSims = 5,
    # gas mole fractions, will be normalized
    # adding CO2 and H2 on their own 
    initialGasMoleFractions = {"H2O": 10,
                               "CO2": 0.05,
                               "H2": 190,
                               "CO": 1,
                               "CH4": 0.8,
                               "N2": 0.85,
                               "O2": 0.15,
                               "NH3": 0.1,},
    # surface coverage, this is proportion, so this is wide open
    initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "H2":0.99,},
    terminationTime=(1000, 's'),
    #terminationRateRatio=0.01,
)
surfaceReactor(
    # range of temps
    temperature=[(650,'K'),(1200,'K')],
    # only works with one pressure, more stuff happening when higher
    initialPressure=(100.0, 'bar'),
    # number of simulations
    nSims = 5,
    # gas mole fractions, will be normalized
    # adding CO2 and H2 on their own 
    initialGasMoleFractions = {"H2O": 0,
                               "CO2": 0,
                               "H2": 8,
                               "CO": 2,
                               "CH4": 0,
                               "N2": 0,
                               "O2": 0,},
    # surface coverage, this is proportion, so this is wide open
    initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "H2":0.99,},
    terminationTime=(1000, 's'),
    #terminationRateRatio=0.01,
)
surfaceReactor(
    # range of temps
    temperature=[(650,'K'),(1200,'K')],
    # only works with one pressure, more stuff happening when higher
    initialPressure=(100.0, 'bar'),
    # number of simulations
    nSims = 5,
    # gas mole fractions, will be normalized
    # adding CO2 and H2 on their own 
    initialGasMoleFractions = {"H2O": 0,
                               "CO2": 2,
                               "H2": 8,
                               "CO": 0,
                               "CH4": 0,
                               "N2": 0,
                               "O2": 0,},
    # surface coverage, this is proportion, so this is wide open
    initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "H2":0.99,},
    terminationTime=(1000, 's'),
    #terminationRateRatio=0.01,
)
surfaceReactor(
    # range of temps
    temperature=[(650,'K'),(1200,'K')],
    # only works with one pressure, more stuff happening when higher
    initialPressure=(100.0, 'bar'),
    # number of simulations
    nSims = 5,
    # gas mole fractions, will be normalized
    # adding CO2 and H2 on their own 
    initialGasMoleFractions = {"H2O": 0,
                               "CO2": 0,
                               "H2": 9,
                               "CO": 0,
                               "CH4": 1,
                               "N2": 0,
                               "O2": 0,},
    # surface coverage, this is proportion, so this is wide open
    initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "H2":0.99,},
    terminationTime=(1000, 's'),
    #terminationRateRatio=0.01,
)
surfaceReactor(
    # range of temps
    temperature=[(650,'K'),(1200,'K')],
    # only works with one pressure, more stuff happening when higher
    initialPressure=(100.0, 'bar'),
    # number of simulations
    nSims = 5,
    # gas mole fractions, will be normalized
    # adding CO2 and H2 on their own 
    initialGasMoleFractions = {"H2O": 0,
                               "CO2": 0,
                               "H2": 9,
                               "CO": 0,
                               "NH3": 1,
                               "N2": 0,
                               "O2": 0,},
    # surface coverage, this is proportion, so this is wide open
    initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "H2":0.99,},
    terminationTime=(1000, 's'),
    #terminationRateRatio=0.01,
)

### Pressure 3

# Reaction systems
surfaceReactor(
    # range of temps
    temperature=[(650,'K'),(1200,'K')],
    # only works with one pressure, more stuff happening when higher
    initialPressure=(0.01, 'bar'),
    # number of simulations
    nSims = 5,
    # gas mole fractions, will be normalized
    # adding CO2 and H2 on their own 
    initialGasMoleFractions = {"H2O": 10,
                               "CO2": 0.05,
                               "H2": 190,
                               "CO": 1,
                               "CH4": 0.8,
                               "N2": 0.85,
                               "O2": 0.15,
                               "NH3": 0.1,},
    # surface coverage, this is proportion, so this is wide open
    initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "H2":0.99,},
    terminationTime=(1000, 's'),
    #terminationRateRatio=0.01,
)
surfaceReactor(
    # range of temps
    temperature=[(650,'K'),(1200,'K')],
    # only works with one pressure, more stuff happening when higher
    initialPressure=(0.01, 'bar'),
    # number of simulations
    nSims = 5,
    # gas mole fractions, will be normalized
    # adding CO2 and H2 on their own 
    initialGasMoleFractions = {"H2O": 0,
                               "CO2": 0,
                               "H2": 8,
                               "CO": 2,
                               "CH4": 0,
                               "N2": 0,
                               "O2": 0,},
    # surface coverage, this is proportion, so this is wide open
    initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "H2":0.99,},
    terminationTime=(1000, 's'),
    #terminationRateRatio=0.01,
)
surfaceReactor(
    # range of temps
    temperature=[(650,'K'),(1200,'K')],
    # only works with one pressure, more stuff happening when higher
    initialPressure=(0.01, 'bar'),
    # number of simulations
    nSims = 5,
    # gas mole fractions, will be normalized
    # adding CO2 and H2 on their own 
    initialGasMoleFractions = {"H2O": 0,
                               "CO2": 2,
                               "H2": 8,
                               "CO": 0,
                               "CH4": 0,
                               "N2": 0,
                               "O2": 0,},
    # surface coverage, this is proportion, so this is wide open
    initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "H2":0.99,},
    terminationTime=(1000, 's'),
    #terminationRateRatio=0.01,
)
surfaceReactor(
    # range of temps
    temperature=[(650,'K'),(1200,'K')],
    # only works with one pressure, more stuff happening when higher
    initialPressure=(0.01, 'bar'),
    # number of simulations
    nSims = 5,
    # gas mole fractions, will be normalized
    # adding CO2 and H2 on their own 
    initialGasMoleFractions = {"H2O": 0,
                               "CO2": 0,
                               "H2": 9,
                               "CO": 0,
                               "CH4": 1,
                               "N2": 0,
                               "O2": 0,},
    # surface coverage, this is proportion, so this is wide open
    initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "H2":0.99,},
    terminationTime=(1000, 's'),
    #terminationRateRatio=0.01,
)
surfaceReactor(
    # range of temps
    temperature=[(650,'K'),(1200,'K')],
    # only works with one pressure, more stuff happening when higher
    initialPressure=(0.01, 'bar'),
    # number of simulations
    nSims = 5,
    # gas mole fractions, will be normalized
    # adding CO2 and H2 on their own 
    initialGasMoleFractions = {"H2O": 0,
                               "CO2": 0,
                               "H2": 9,
                               "CO": 0,
                               "NH3": 1,
                               "N2": 0,
                               "O2": 0,},
    # surface coverage, this is proportion, so this is wide open
    initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "H2":0.99,},
    terminationTime=(1000, 's'),
    #terminationRateRatio=0.01,
)


simulator(
    atol=1e-18,
    rtol=1e-12,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.5,
    toleranceInterruptSimulation=0.5,
    #toleranceBranchReactionToCore=0.001,
    #branchingIndex=0.5,
    #branchingRatioMax=1.0,
    #maximumEdgeSpecies=500000,
)

options(
    units='si',
    generateOutputHTML=True,
    generatePlots=False, # Enable to make plots of core and edge size etc. But takes a lot of the total runtime!
    saveEdgeSpecies=False,
    saveSimulationProfiles=False,
    generateCanteraYAML2={'saveInterval': 1, 'saveEdge': False},
)

generatedSpeciesConstraints(
    allowed=['input species', 'seed mechanisms', 'reaction libraries'],
    maximumCarbonAtoms=2,
    maximumOxygenAtoms=2,
    maximumNitrogenAtoms=2,
    maximumSurfaceSites=2,
    allowSingletO2=True,
)








