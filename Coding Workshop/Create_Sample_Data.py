import random
import pandas as pd
from datetime import datetime, timedelta
from tqdm import tqdm

random.seed(42)

# --------------------------------------------- Note ---------------------------------------------
# Resistance and pathogen data for the outbreak was taken from: 
# https://pubchem.ncbi.nlm.nih.gov/protein/P0A5N1#section=Enzyme-Commission-%28EC%29-Number 
# https://www.ncbi.nlm.nih.gov/pathogens/refgene/#
# ------------------------------------------------------------------------------------------------

# Configuration
NUM_PATIENTS = 800
OUTBREAK_CASES = 700
START_DATE = datetime(2025, 1, 1)
DURATION = 200
NON_RESISTANT_RATE = 0.2 # 20% of background cases will be non-resistant 

# --------------- PATHOGEN DATA ---------------
BACKGROUND_PATHOGENS = [
    'Escherichia coli',
    'Klebsiella pneumoniae',
    'Pseudomonas aeruginosa',
    'Methicillin-Resistant Staphylococcus aureus'
]

OUTBREAK_PATHOGEN = 'Mycobacterium tuberculosis var. bovis AF2122/97'

# --------------- RESISTANCE GENE PROFILES ---------------
RESISTANCE_GENES = {
    'Mycobacterium tuberculosis var. bovis AF2122/97': ['aac(2)-Ic'],
    'Escherichia coli': ['blaCTX-M', 'aac(3)-IIa'],
    'Klebsiella pneumoniae': ['blaKPC', 'blaOXA-48'],
    'Pseudomonas aeruginosa': ['mexA', 'oprD'],
    'Methicillin-Resistant Staphylococcus aureus': ['mecA', 'blaZ']
}

# --------------- CLINICAL METADATA ---------------
SEX = ['Male', 'Female']
OUTCOMES = ['Recovered', 'Deceased', 'Ongoing']
UNITS = ['ICU-1', 'ICU-2', 'ICU-3', 'ED', 'Med-Surg', 'CCU', 'PCU']
SAMPLE_SITES = ['Sputum', 'Lung Biopsy', 'Blood', 'CSF']
TREATMENTS = ['Isoniazid', 'Rifampicin', 'Ethambutol', 'Streptomycin', 'Kanamycin', 'None']

# --------------- MICROBIAL / GENETIC DATA ---------------
SEQ_PLATFORMS = ['Illumina MiSeq', 'Oxford Nanopore', 'PacBio']
MEDIUMS = ['7H9 broth', 'Lowenstein-Jensen', 'Middlebrook agar']
PROTOCOLS = ['Nextera XT', 'TruSeq DNA', 'Ligation Sequencing Kit']
GENE_CLUSTERS = {
    'ESBL + Carbapenemase': ['blaCTX-M', 'blaKPC'],
    'MRSA Cluster': ['mecA', 'blaZ'],
    'Pseudomonas Cluster': ['mexA', 'oprD']
}

STRAIN_PREFIX = {
    'Escherichia coli': 'ECOLI',
    'Klebsiella pneumoniae': 'KPN',
    'Pseudomonas aeruginosa': 'PAE',
    'Methicillin-Resistant Staphylococcus aureus': 'MRSA'
}

Strain_Counters = {
    'ECOLI': 1,
    'KPN': 1,
    'PAE': 1,
    'MRSA': 1
}

# ---------------------------------- OUTBREAK TIMING -----------------------------------
# Simulate two-wave outbreak: peaks around day 45 and day 120
Day_Weights = []
for day in range(DURATION):
    Wave1 = max(0, 30 - abs(day - 45))  # Gradual rise/fall over 30 days before/after
    Wave2 = max(0, 30 - abs(day - 120))
    Total_Weight = Wave1 + Wave2
    Day_Weights.append(Total_Weight)

# Generate outbreak collection days for each outbreak case
Outbreak_Collection_Days = random.choices(
    population = range(DURATION),
    weights = Day_Weights,
    k = OUTBREAK_CASES
)
# ---------------------------------- DATA GENERATION -----------------------------------

print('Processing Clinical and Laboratory Dataset .........')

# Data generation
Clinical_Data = []
Laboratory_Data = [] 

# Randomly select indices for outbreak cases
Random_Indicies = set(random.sample(range(NUM_PATIENTS), OUTBREAK_CASES))
Outbreak_Day_Iterator = iter(Outbreak_Collection_Days)

for i in tqdm(range(NUM_PATIENTS)):
    Patient_ID = 2200 + i
    Sex = random.choice(SEX)
    Age = random.randint(18, 85)
    Unit = random.choice(UNITS)
    Sample_Site = random.choice(SAMPLE_SITES)
    Duration_Stay = random.randint(5, 60)
    Treatment = random.choice(TREATMENTS)
    Platform = random.choice(SEQ_PLATFORMS)
    Growth_Medium = random.choice(MEDIUMS)
    Libary_Prep = random.choice(PROTOCOLS)
    Incubation_Time = random.randint(18, 48)
    Genome_Coverage = round(random.uniform(90.0, 99.9), 2)
    Quality_Score = round(random.uniform(30.0, 40.0), 2)

    if i in Random_Indicies:
        # Outbreak case
        Strain = OUTBREAK_PATHOGEN
        Diagnosis = 'Drug-resistant TB'
        Genes = RESISTANCE_GENES[OUTBREAK_PATHOGEN]
        MLST = 'ST42'
        Outbreak_Day = next(Outbreak_Day_Iterator)
        Collection_Date = START_DATE + timedelta(days=Outbreak_Day)

        Expression_Levels = {gene: round(random.uniform(450.0, 1000.0), 2) for gene in Genes}
        Strain_ID = 'bTB-R1'

    else:
        # Background infection
        Strain = random.choice(BACKGROUND_PATHOGENS)
        # For background infections over a wide time window (before + during)
        Offset = random.randint(0, DURATION - 1)
        Collection_Date = START_DATE + timedelta(days = Offset)

        Diagnosis = f'Infection - {Strain}'
        # Decide if its a non-resistant strain 
        Non_Resistant = random.random() < NON_RESISTANT_RATE
        Prefix = STRAIN_PREFIX[Strain]

        if Non_Resistant: 
            Genes = [] 
            Expression_Levels = {} 

            Strain_ID = f"{Prefix}-{Strain_Counters[Prefix]}"
        else: 
            # 30% chance of cluster genes 
            Cluster_Genes = [] 
            if random.random() < 0.3:
                if Strain == 'Escherichia coli' or Strain == 'Klebsiella pneumoniae':
                    Cluster_Genes = GENE_CLUSTERS['ESBL + Carbapenemase']
                elif Strain == 'Methicillin-Resistant Staphylococcus aureus':
                    Cluster_Genes = GENE_CLUSTERS['MRSA Cluster']
                elif Strain == 'Pseudomonas aeruginosa':
                    Cluster_Genes = GENE_CLUSTERS['Pseudomonas Cluster']
            
            Strain_ID = f"{Prefix}-R{Strain_Counters[Prefix]}"
        
        ### Maybe you cna change it if it has the same resistance genes as another test then label as same strain
        Strain_Counters[Prefix] += 1
        
        Available_Genes = RESISTANCE_GENES[Strain]
        if Cluster_Genes:
            Genes = Cluster_Genes
        else:
            Num_Genes = 0
            Genes = []
        
        Expression_Levels = {gene: round(random.uniform(10.0, 300.00), 2) for gene in Genes}

    if Num_Genes == 0: 
        Avg_Expression = 0
    else: 
        Avg_Expression = sum(Expression_Levels.values()) / len(Expression_Levels)
        
    if Avg_Expression > 700:
        Outcome = random.choices(OUTCOMES, weights=[0.2, 0.7, 0.1])[0]
    elif Avg_Expression > 500:
        Outcome = random.choices(OUTCOMES, weights=[0.10, 0.45, 0.45])[0]
    else:
        Outcome = random.choices(OUTCOMES, weights=[0.7, 0.1, 0.2])[0]


    MLST = f'ST{random.randint(10, 50)}'


    Clinical_Data.append({ 
        'Patient ID': Patient_ID, 
        'Sex': Sex, 
        'Age': Age, 
        'Unit': Unit, 
        'Diagnosis': Diagnosis, 
        'Treatment': Treatment, 
        'Outcome': Outcome, 
        'Hospital_Duration': Duration_Stay
    })

    Laboratory_Data.append({
        'Patient ID': Patient_ID, 
        'Collection Date': Collection_Date, 
        'Sample Site': Sample_Site, 
        'Strain': Strain, 
        'Strain ID': Strain_ID,
        'Resistance Genes': Genes, 
        'Gene Expression Levels': Expression_Levels,
        'MLST': MLST, 
        'Genome Coverage': Genome_Coverage, 
        'Sample Quality Score': Quality_Score, 
        'Sequencing Platform': Platform, 
        'Growth Media': Growth_Medium, 
        'Library Prep Protocol': Libary_Prep, 
        'Incubation Time [HR]': Incubation_Time
    })

# Convert to DataFrame and save
Clinical_DF = pd.DataFrame(Clinical_Data)
Clinical_DF.to_csv('Clinical_Data.csv', index=False)
print('Clinical Dataset saved as "Clinical_Dataset.csv"')

Lab_DF = pd.DataFrame(Laboratory_Data)
Lab_DF.to_csv('Lab_Processing_Data.csv', index = False)
print('Laboratory Dataset saved as "Lab_Processing_Data.csv"')

# Combined resistance protein sequences (FASTA format)
print('Processing Resistance Data .........')
Resistance_Data = {
    "aac(2)-Ic": (
        ">aac(2)-Ic | Mycobacterium tuberculosis variant bovis AF2122/97\n"
        "MHTQVHTARLVHTADLDSETRQDIRQMVTGAFAGDFTETDWEHTLGGMHALIWHHGAIIAHAAVIQRRLIYRGNALRCGYVEGVAVRADWRGQRL"
        "VSALLDAVEQVMRGAYQLGALSSSARARRLYASRGWLPWHGPTSVLAPTGPVRTPDDDGTVFVLPIDISLDTSAELMCDWRAGDVW"
    ),
    "blaCTX-M": (
        ">blaCTX-M | E. coli extended-spectrum beta-lactamase\n"
        "MKKTAIAIAVALAGFATVAQAAPTASVVVVTGAAQSRVTAATSTLQKLLGRIVTGYVNNFNDPTTTGTPVGIVQ"
        "TQHSTKPLRLLTSDTRQQLDLAGQKPLPSGSLRVEISFALAGVDQELTTTQKAPLS"
    ),
    "aac(3)-IIa": (
        ">aac(3)-IIa | E. coli aminoglycoside acetyltransferase\n"
        "MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWMTDNMLRAAMRFGELFQHDIVQDAAVLKKCGAWLS"
    ),
    "blaKPC": (
        ">blaKPC | Klebsiella pneumoniae carbapenemase\n"
        "MKKTAIAIAVALAGFATVAQAAPTASVVVVTGAAQSRVTAATSTLQKLLGRIVTGYVNNFNDPTTTGTPV"
    ),
    "blaOXA-48": (
        ">blaOXA-48 | Klebsiella pneumoniae oxacillinase\n"
        "MKIKTGARLLFTAAALAVALAPAIAAQITTAQQLLDATTHVDFVPSLGDPAALRRLINPQGLLPKPETDF"
    ),
    "mexA": (
        ">mexA | Pseudomonas aeruginosa multidrug efflux protein\n"
        "MNKTLLSLALLLAVSAVAGQAAFAEQDAGADLISLTDTPETVAEAFNLLAGHPVLHGRLLQNSDGSVF"
    ),
    "oprD": (
        ">oprD | Pseudomonas aeruginosa outer membrane porin\n"
        "MKTVLFAGLVATAGGSAQQGQVSFTRGEKVYQPTGQTSFNSISADYQKVLDSDGTVVRANQ"
    ),
    "mecA": (
        ">mecA | MRSA penicillin-binding protein PBP2a\n"
        "MLLKKTKIKIVPLILVALIVVLVLIGTTATAIEGNLVTPIWAAITILGAKVVNSPGGYP"
    ),
    "blaZ": (
        ">blaZ | MRSA beta-lactamase\n"
        "MKNTLSFASKVFGAACATYTATKKLKGKTLTNSDQVRDIKKLIEKQGADTVEAAAKEAA"
    )
}

# Save to a .txt file with FASTA format content
with open("RESISTANCE_GENES.txt", "w") as f:
    for Seq in tqdm(Resistance_Data.values()):
        f.write(Seq + "\n")

print('Resistance gene sequences saved to "RESISTANCE_GENES.txt"')






print("Processing Cattle Farm Data ........")

Cattle_Test_Records = []
FARMS = ['FARM_001', 'FARM_002', 'FARM_003', 'FARM_004', 'FARM_005']
RESISTANT_TB_FARMS = ['FARM_001', 'FARM_004']
HERDS_PER_FARM = 10
CATTLE_PER_HERD = 10
BREEDS = ['Holstein', 'Jersey', 'Angus', 'Guernsey', 'Brown Swiss']
TEST_TYPES = ['PPD Skin Test', 'Gamma Interferon Test']

for Farm_ID in tqdm(FARMS):
    for Herd_Num in range(1, HERDS_PER_FARM + 1):
        Herd_ID = f'H{Herd_Num:02}'
        for Cow_Num in range(1, CATTLE_PER_HERD + 1):
            Cow_ID = f'C{Cow_Num:03}'
            Age_Months = random.randint(6, 120)
            Cow_Sex = random.choice(['Male', 'Female'])
            Breed = random.choice(BREEDS)
            Body_Score = round(random.uniform(1.0, 5.0), 1)
            Milk_Production = round(random.uniform(10.0, 35.0), 1) if Cow_Sex == 'Female' else None
            Test_Type = random.choice(TEST_TYPES)
            Test_Date = START_DATE - timedelta(days=random.randint(0, 60))

            Gene_Tested = "aac(2)-Ic"
            aac_Present = False
            Expression = None
            Cow_Strain_ID = None
            Zoonotic_Risk = 'No'

            if Farm_ID in RESISTANT_TB_FARMS:
                # TB-positive? (45% chance)
                TB_Postive = random.random() < 0.45

                # TB-positive logic
                if TB_Postive:

                    # 60% get resistant gene
                    if random.random() < 0.6:
                        aac_Present = True
                        Expression = round(random.uniform(100.0, 1000.0), 2)
                        Cow_Strain_ID = 'bTB-R1'
                        Zoonotic_Risk = 'Yes'
                    else:
                        Cow_Strain_ID = 'bTB-WT'
                else:
                    Cow_Strain_ID = 'bTB-WT'
            else:
                # TB-positive? (5% chance)
                TB_Postive = random.random() < 0.05
                # Possible false-positive detection in healthy animals
                if random.random() < 0.01:
                    aac_Present = True
                    Expression = round(random.uniform(10.0, 300.0), 2)
                    Cow_Strain_ID = 'bTB-R1'
                    Zoonotic_Risk = 'Yes'

            # Diagnosis
            if TB_Postive and aac_Present:
                Diagnosis = 'Drug-resistant bTB'
                Test_Result = 'Positive'
            elif TB_Postive:
                Diagnosis = 'bTB'
                Test_Result = 'Positive'
            elif aac_Present:
                Diagnosis = 'Carrier / Monitored'
                Test_Result = 'Suspicious'
            else:
                Diagnosis = 'Healthy'
                Test_Result = 'Negative'

            Follow_Up = Test_Result in ['Positive', 'Suspicious']
            Notes = 'Isolated' if Diagnosis == 'Drug-resistant bTB' else (
                'Retest Scheduled' if Follow_Up else 'Healthy')

            Cattle_Test_Records.append({
                'Farm ID': Farm_ID,
                'Herd ID': Herd_ID,
                'Cow ID': Cow_ID,
                'Age (Months)': Age_Months,
                'Sex': Cow_Sex,
                'Breed': Breed,
                'Body Condition Score': Body_Score,
                'Milk Production (L/day)': Milk_Production,
                'Test Date': Test_Date.strftime('%Y-%m-%d'),
                'Test Type': Test_Type,
                'Test Result': Test_Result,
                'Diagnosis': Diagnosis,
                'Strain ID': Cow_Strain_ID,
                'Gene Target Tested': Gene_Tested,
                'aac(2)-Ic Present': aac_Present,
                'aac(2)-Ic Expression': Expression,
                'Follow-up Needed': Follow_Up,
                'Notes': Notes
                })

# Save new cattle dataset
Cattle_Farm_DF = pd.DataFrame(Cattle_Test_Records)
Cattle_Farm_DF.to_csv("Cattle_Farm_Data.csv", index=False)
print('Cattle Farm TB Dataset saved as "Cattle_Farm_Data.csv"')
