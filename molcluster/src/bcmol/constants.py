"""Structure constants for the whole project were provided here."""

ATOMIC_NUMBERS = {
    'H': 1, 'D': 1, 'T': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20,
    'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
    'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40,
    'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
    'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60,
    'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70,
    'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80,
    'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90,
    'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100,
    'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109,
    'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118
}

ELEMENTS = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
    21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
    41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
    61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
    71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
    91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
    101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds',
    111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'
}

# GENERAL_ATOMIC_NUMBERS = {
#     'H': 1,
#     'C': 6,
#     'N': 7,
#     'O': 8,
#     'F': 9,
#     'Si': 14,  # Supported element in MMFF94 and Dock6 but not in Vina
#     'P': 15,
#     'S': 16,
#     'Cl': 17,
#     'Se': 34,  # Supported element in Dock6 but not in MMFF94 or Vina
#     'Br': 35,
#     'I': 53,
# }

GENERAL_ELEMENTS = ['H', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Se', 'Br', 'I']
MMFF94_ELEMENTS = ['H', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Br', 'I', 'Fe', 'Cu', 'Zn', 'Mg', 'Na', 'K', 'Ca']
METAL_ELEMENTS = ['Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
                  'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',
                  'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Cs', 'Ba', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os',
                  'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'Fr', 'Ra', 'Lr', 'Ho']
VINA_ELEMENTS = ['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I'] + METAL_ELEMENTS
DOCK6_ELEMENTS = ['H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Se', 'Br', 'I'] + METAL_ELEMENTS
GENERAL_ATOMIC_NUMBERS = [ATOMIC_NUMBERS[k] for k in GENERAL_ELEMENTS]
MMFF94_ATOMIC_NUMBERS = [ATOMIC_NUMBERS[k] for k in MMFF94_ELEMENTS]
METAL_ATOMIC_NUMBERS = [ATOMIC_NUMBERS[k] for k in METAL_ELEMENTS]
VINA_ATOMIC_NUMBERS = [ATOMIC_NUMBERS[k] for k in VINA_ELEMENTS]
DOCK6_ATOMIC_NUMBERS = [ATOMIC_NUMBERS[k] for k in DOCK6_ELEMENTS]

AA = 'ACDEFGHIKLMNPQRSTVWY'  # IUPAC.protein.letters
BASES = ['C', 'T', 'G', 'U', 'A', 'I']
ALL_BASES = BASES + ['%s5' % _b for _b in BASES] + ['%s3' % _b for _b in BASES]
DBASES = ['D%s' % _b for _b in ALL_BASES]
RBASES = ['R%s' % _b for _b in ALL_BASES] + BASES

# Add terminal base names (sometimes used by forcefields)
for _nucleic_list in DBASES, RBASES:
    for _base in list(_nucleic_list):
        _nucleic_list.append(_base + '3')
        _nucleic_list.append(_base + '5')

AA321 = dict(ALA='A', ASX='B', CYS='C', ASP='D', GLU='E', PHE='F', GLY='G', HIS='H', ILE='I',
             XLE='J', LYS='K', LEU='L', MET='M', ASN='N', PYL='O', PRO='P', GLN='Q', ARG='R',
             SER='S', THR='T', SEC='U', VAL='V', TRP='W', XAA="X", TYR='Y', GLX="Z")
AA123 = {v: k for k, v in AA321.items()}

# AMINO_NAMES = {
#     "ALA": "Alanine",
#     "ARG": "Arginine",
#     "ASH": "Neutral ASP",
#     "ASN": "Asparagine",
#     "ASP": "Aspartic acid",
#     "ASX": "ASP/ASN ambiguous",
#     "CYS": "Cysteine",
#     "CYX": "Cystine, SS-bonded CYS",
#     "CYM": "Cysteine anion",
#     "GLH": "Neutral GLU",
#     "GLN": "Glutamine",
#     "GLU": "Glutamic acid",
#     "GLX": "GLU/GLN ambiguous",
#     "GLY": "Glycine",
#     "HIS": "Histidine",
#     "HIE": "Histidine epsilon tautomer; Neutral HIS, proton HE2 present",
#     "HID": "Histidine delta tautomer; Neutral HIS, proton HD1 present",
#     "HIP": "Histidine ion",
#     "HSD": "",
#     "HSE": "L-Homoserine",
#     "HSP": "",
#     "ILE": "Isoleucine",
#     "LEU": "Leucine",
#     "LYN": "Neutral LYS",
#     "LYS": "Lysine",
#     "MET": "Methionine",
#     "PHE": "Phenylalanine",
#     "PRO": "Proline",
#     "SER": "Serine",
#     "THR": "Threonine",
#     "TRP": "Tryptophan",
#     "TYM": "Negative TYR",
#     "TYR": "Tyrosine",
#     "VAL": "Valine",
#     # "UNK": "Undetermined",
#     # "SEP": "Phosphorylated Serine",
#     # "MEP": "Phosphorylated Methionine"  # TODO: add more non-standard residues
# }

AMINO_NAMES = AA321

BACKBONES = {'dna': set(("P OP1 OP2 O5' O4' C5' C4' C3' O3' C2' C1' H1' H2'' H2' H3' H4' H5' H5'' "
                         "HO5' HO3'").split()),
             'protein': set("N CA C O OXT H HA HA2 HA3 H2 H3".split())}

CAPS = {'ACE': 'Acetyl',
        'NME': 'N-methyl amide'}

WATER_NAMES = [
    'SOL', 'WAT', 'HOH', 'H2O', 'W', 'DOD', 'D3O', 'TIP3', 'TIP4', 'SPC',
]

# all chemical components with the word "ion" in their name, Sep 2016
#
# SET SESSION group_concat_max_len = 1000000;
# SELECT GROUP_CONCAT(id_ ORDER BY id_ ASC SEPARATOR '", "') from
# (
#     SELECT count(obj_id) as c, id_
#     FROM pdb.chem_comp WHERE name LIKE "% ION%"
#     GROUP BY id_
# ) AS t1;
ION_NAMES = [
    '118', '119', '1AL', '1CU', '2FK', '2HP', '2OF', '3CO',
    '3MT', '3NI', '3OF', '3P8', '4MO', '4PU', '543', '6MO', 'ACT', 'AG', 'AL',
    'ALF', 'AM', 'ATH', 'AU', 'AU3', 'AUC', 'AZI', 'BA', 'BCT', 'BEF', 'BF4', 'BO4',
    'BR', 'BS3', 'BSY', 'CA', 'CAC', 'CD', 'CD1', 'CD3', 'CD5', 'CE', 'CHT', 'CL',
    'CO', 'CO3', 'CO5', 'CON', 'CR', 'CS', 'CSB', 'CU', 'CU1', 'CU3', 'CUA', 'CUZ',
    'CYN', 'DME', 'DMI', 'DSC', 'DTI', 'DY', 'E4N', 'EDR', 'EMC', 'ER3', 'EU',
    'EU3', 'F', 'FE', 'FE2', 'FPO', 'GA', 'GD3', 'GEP', 'HAI', 'HG', 'HGC', 'IN',
    'IOD', 'IR', 'IR3', 'IRI', 'IUM', 'K', 'KO4', 'LA', 'LCO', 'LCP', 'LI', 'LU',
    'MAC', 'MG', 'MH2', 'MH3', 'MLI', 'MLT', 'MMC', 'MN', 'MN3', 'MN5', 'MN6',
    'MO1', 'MO2', 'MO3', 'MO4', 'MO5', 'MO6', 'MOO', 'MOS', 'MOW', 'MW1', 'MW2',
    'MW3', 'NA', 'NA2', 'NA5', 'NA6', 'NAO', 'NAW', 'NCO', 'NET', 'NH4', 'NI',
    'NI1', 'NI2', 'NI3', 'NO2', 'NO3', 'NRU', 'O4M', 'OAA', 'OC1', 'OC2', 'OC3',
    'OC4', 'OC5', 'OC6', 'OC7', 'OC8', 'OCL', 'OCM', 'OCN', 'OCO', 'OF1', 'OF2',
    'OF3', 'OH', 'OS', 'OS4', 'OXL', 'PB', 'PBM', 'PD', 'PDV', 'PER', 'PI', 'PO3',
    'PO4', 'PR', 'PT', 'PT4', 'PTN', 'RB', 'RH3', 'RHD', 'RU', 'SB', 'SCN', 'SE4',
    'SEK', 'SM', 'SMO', 'SO3', 'SO4', 'SR', 'T1A', 'TB', 'TBA', 'TCN', 'TEA', 'TH',
    'THE', 'TL', 'TMA', 'TRA', 'UNX', 'V', 'VN3', 'VO4', 'W', 'WO5', 'Y1', 'YB',
    'YB2', 'YH', 'YT3', 'ZCM', 'ZN', 'ZN2', 'ZN3', 'ZNO', 'ZO3',
    # additional ion names
    'OHX', 'Na+', 'Cl-'
]

# all chemical components with the word "%saccharide%" in their type, Sep 2016
#
# SET SESSION group_concat_max_len = 1000000;
# select GROUP_CONCAT(id_ ORDER BY id_ ASC SEPARATOR '", "') from
# (
#     SELECT count(obj_id), id_
#     FROM pdb.chem_comp WHERE type like "%SACCHARIDE%"
#     GROUP BY id_
# ) AS t1;
SACCHARIDE_NAMES = [
    '045', '0AT', '0BD', '0MK', '0NZ', '0TS', '0V4', '0XY', '0YT', '10M',
    '147', '149', '14T', '15L', '16G', '18T', '18Y', '1AR', '1BW', '1GL', '1GN',
    '1JB', '1LL', '1NA', '1S3', '26M', '26Q', '26R', '26V', '26W', '26Y', '27C',
    '289', '291', '293', '2DG', '2F8', '2FG', '2FL', '2FP', '2GL', '2M4', '2M5',
    '32O', '34V', '3CM', '3DO', '3DY', '3FM', '3LR', '3MF', '3MG', '3SA', '3ZW',
    '46D', '46M', '46Z', '48Z', '4CQ', '4GC', '4NN', '50A', '5DI', '5GF', '5MM',
    '5RP', '5SA', '5SP', '64K', '6PG', '6SA', '7JZ', '7SA', 'A1Q', 'A2G', 'AAB',
    'AAL', 'AAO', 'ABC', 'ABD', 'ABE', 'ABF', 'ABL', 'ACG', 'ACI', 'ACR', 'ACX',
    'ADA', 'ADG', 'ADR', 'AF1', 'AFD', 'AFL', 'AFO', 'AFP', 'AFR', 'AGC', 'AGH',
    'AGL', 'AHR', 'AIG', 'ALL', 'ALX', 'AMU', 'AOG', 'AOS', 'ARA', 'ARB', 'ARE',
    'ARI', 'ASG', 'ASO', 'AXP', 'AXR', 'B0D', 'B16', 'B2G', 'B4G', 'B6D', 'B8D',
    'B9D', 'BBK', 'BCD', 'BDG', 'BDP', 'BDR', 'BEM', 'BFP', 'BGC', 'BGL', 'BGP',
    'BGS', 'BHG', 'BMA', 'BMX', 'BNG', 'BNX', 'BOG', 'BRI', 'BXF', 'BXP', 'BXX',
    'BXY', 'C3X', 'C4X', 'C5X', 'CAP', 'CBI', 'CBK', 'CBS', 'CDR', 'CEG', 'CGF',
    'CHO', 'CR1', 'CR6', 'CRA', 'CT3', 'CTO', 'CTR', 'CTT', 'D6G', 'DAF', 'DAG',
    'DDA', 'DDB', 'DDL', 'DEL', 'DFR', 'DFX', 'DG0', 'DGC', 'DGD', 'DGM', 'DGS',
    'DIG', 'DLF', 'DLG', 'DMU', 'DNO', 'DOM', 'DP5', 'DQQ', 'DQR', 'DR2', 'DR3',
    'DR4', 'DRI', 'DSR', 'DT6', 'DVC', 'E4P', 'E5G', 'EAG', 'EBG', 'EBQ', 'EGA',
    'EJT', 'EPG', 'ERE', 'ERI', 'F1P', 'F1X', 'F6P', 'FBP', 'FCA', 'FCB', 'FCT',
    'FDP', 'FDQ', 'FFC', 'FIX', 'FMO', 'FRU', 'FSI', 'FU4', 'FUB', 'FUC', 'FUD',
    'FUL', 'FXP', 'G16', 'G1P', 'G2F', 'G3I', 'G4D', 'G4S', 'G6D', 'G6P', 'G6S',
    'GAC', 'GAD', 'GAL', 'GC1', 'GC4', 'GCD', 'GCN', 'GCO', 'GCS', 'GCT', 'GCU',
    'GCV', 'GCW', 'GCX', 'GE1', 'GFG', 'GFP', 'GIV', 'GL0', 'GL2', 'GL5', 'GL6',
    'GL7', 'GL9', 'GLA', 'GLB', 'GLC', 'GLD', 'GLF', 'GLG', 'GLO', 'GLP', 'GLS',
    'GLT', 'GLW', 'GMH', 'GN1', 'GNX', 'GP1', 'GP4', 'GPH', 'GPM', 'GQ1', 'GQ2',
    'GQ4', 'GS1', 'GS4', 'GSA', 'GSD', 'GTE', 'GTH', 'GTK', 'GTR', 'GTZ', 'GU0',
    'GU1', 'GU2', 'GU3', 'GU4', 'GU5', 'GU6', 'GU8', 'GU9', 'GUF', 'GUP', 'GUZ',
    'GYP', 'GYV', 'H2P', 'HDL', 'HMS', 'HS2', 'HSD', 'HSG', 'HSH', 'HSJ', 'HSQ',
    'HSR', 'HSU', 'HSX', 'HSY', 'HSZ', 'IAB', 'IDG', 'IDR', 'IDS', 'IDT', 'IDU',
    'IDX', 'IDY', 'IMK', 'IN1', 'IPT', 'ISL', 'KBG', 'KD2', 'KDA', 'KDM', 'KDO',
    'KFN', 'KO1', 'KO2', 'KTU', 'L6S', 'LAG', 'LAI', 'LAK', 'LAO', 'LAT', 'LB2',
    'LBT', 'LCN', 'LDY', 'LGC', 'LGU', 'LM2', 'LMT', 'LMU', 'LOG', 'LOX', 'LPK',
    'LSM', 'LTM', 'LVZ', 'LXB', 'LXZ', 'M1F', 'M3M', 'M6P', 'M8C', 'MA1', 'MA2',
    'MA3', 'MAB', 'MAG', 'MAL', 'MAN', 'MAT', 'MAV', 'MAW', 'MBG', 'MCU', 'MDA',
    'MDM', 'MDP', 'MFA', 'MFB', 'MFU', 'MG5', 'MGA', 'MGL', 'MLB', 'MMA', 'MMN',
    'MN0', 'MRP', 'MTT', 'MUG', 'MVP', 'MXY', 'N1L', 'N9S', 'NAA', 'NAG', 'NBG',
    'NDG', 'NED', 'NG1', 'NG6', 'NGA', 'NGB', 'NGC', 'NGE', 'NGF', 'NGL', 'NGS',
    'NGY', 'NHF', 'NM6', 'NM9', 'NTF', 'NTO', 'NTP', 'NXD', 'NYT', 'OPG', 'OPM',
    'ORP', 'OX2', 'P3M', 'P53', 'P6P', 'PA5', 'PNA', 'PNG', 'PNW', 'PRP', 'PSJ',
    'PSV', 'PTQ', 'QDK', 'QPS', 'QV4', 'R1P', 'R1X', 'R2B', 'R5P', 'RAA', 'RAE',
    'RAF', 'RAM', 'RAO', 'RAT', 'RB5', 'RBL', 'RCD', 'RDP', 'REL', 'RER', 'RF5',
    'RG1', 'RGG', 'RHA', 'RIB', 'RIP', 'RNS', 'RNT', 'ROB', 'ROR', 'RPA', 'RST',
    'RUB', 'RUU', 'RZM', 'S6P', 'S7P', 'SA0', 'SCR', 'SDD', 'SF6', 'SF9', 'SG4',
    'SG5', 'SG6', 'SG7', 'SGA', 'SGC', 'SGD', 'SGN', 'SGS', 'SHB', 'SHG', 'SI3',
    'SIO', 'SOE', 'SOL', 'SSG', 'SUC', 'SUP', 'SUS', 'T6P', 'T6T', 'TAG', 'TCB',
    'TDG', 'TGK', 'TGY', 'TH1', 'TIA', 'TM5', 'TM6', 'TM9', 'TMR', 'TMX', 'TOA',
    'TOC', 'TRE', 'TYV', 'UCD', 'UDC', 'VG1', 'X0X', 'X1X', 'X2F', 'X4S', 'X5S',
    'X6X', 'XBP', 'XDN', 'XDP', 'XIF', 'XIM', 'XLF', 'XLS', 'XMM', 'XUL', 'XXR',
    'XYP', 'XYS', 'YO5', 'Z3Q', 'Z6J', 'Z9M', 'ZDC', 'ZDM'
]


RESTYPES = dict(
    Protein=set(AMINO_NAMES),
    DNA=set(DBASES),
    RNA=set(RBASES),
    Saccharide=set(SACCHARIDE_NAMES),
    Ion=set(ION_NAMES),
    Water=set(WATER_NAMES),
    Hetero=set(CAPS)
)


def _make_residue_type_dict():
    rdt = {None: 'placeholder'}
    for typename, namelist in RESTYPES.items():
        for resname in namelist:
            rdt[resname] = typename
    return rdt


RESIDUE_TYPES = _make_residue_type_dict()
