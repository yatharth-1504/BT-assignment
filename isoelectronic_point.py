# This is the sequence of the chain A it has been extracted from chain A
sequence = ['MET', 'GLU', 'LEU', 'LYS', 'HIS', 'SER', 'ILE', 'SER', 'ASP', 'TYR', 'THR', 'GLU', 'ALA', 'GLU', 'PHE', 'LEU', 'GLN', 'LEU', 'VAL', 'THR', 'THR', 'ILE', 'CYS', 'ASN', 'ALA', 'ASP', 'THR', 'SER', 'SER', 'GLU', 'GLU', 'GLU', 'LEU', 'VAL', 'LYS', 'LEU', 'VAL', 'THR', 'HIS', 'PHE', 'GLU',
            'GLU', 'MET', 'THR', 'GLU', 'HIS', 'PRO', 'SER', 'GLY', 'SER', 'ASP', 'LEU', 'PRO', 'LYS', 'GLU', 'GLY', 'ASP', 'ASP', 'ASP', 'SER', 'PRO', 'SER', 'GLY', 'ILE', 'VAL', 'ASN', 'THR', 'VAL', 'LYS', 'GLN', 'TRP', 'ARG', 'ALA', 'ALA', 'ASN', 'GLY', 'LYS', 'SER', 'GLY', 'PHE', 'LYS', 'GLN', 'GLY']


# Here we are calculating frequency of a charged residue in the given chain sequence.
def count_of_charged_residue(name):
    cnt = 0
    for i in range(len(sequence)):
        if sequence[i] == name:
            cnt += 1
    return cnt


# Data represents list of dicts. these dicts contain name frequency and pka values of the charged residues.
data = [{"name": "ASP", "count": count_of_charged_residue(
    "ASP"), "pka": 3.65}, {"name": "GLU", "count": count_of_charged_residue("GLU"), "pka": 4.25},
    {"name": "HIS", "count": count_of_charged_residue("HIS"), "pka": 6.00},
    {"name": "LYS", "count": count_of_charged_residue("LYS"), "pka": 10.53},
    {"name": "ARG", "count": count_of_charged_residue("ARG"), "pka": 12.48}]


# These are the critical ph points, they represents ph values where we may experience charge change.
# 2.34 represents pka1 value of GLY, it is been considered since it ends our sequence.
# 9.21 represents pka2 value of MET, it is been considered since it starts our sequence.
critical_ph_points = [2, 3, 4, 5,  8, 11, 13]


# Function to return iso-electronic point of our protien.
# We are iterating over critical_ph_points list and checking at what ph the sign of the net charge changes.
# Then we calculate the wieghted mean of the pka3 value which surround the critical ph obtained.
def cal_iso_pt():
    charge = 0
    for ph in critical_ph_points:
        prev_charge = charge
        charge = {ph < 2.34: +1, 2.34 < ph < 9.21: 0}.get(True, -1)
        for item in data:
            if item["name"] in ["ASP", "GLU"]:
                if ph > item["pka"]:
                    charge -= item["count"]
            if item["name"] in ["ARG", "LYS", "HIS"]:
                if ph < item["pka"]:
                    charge += item["count"]
        if charge*prev_charge <= 0:
            if 3.65 < ph < 4.25:
                return ((3.65*data[0]["count"] + 4.25*data[1]["count"]) / (data[0]["count"] + data[1]["count"]))
            elif 4.25 < ph < 6:
                return ((4.25*data[1]["count"] + 6*data[2]["count"]) / (data[1]["count"] + data[2]["count"]))
            elif 6 < ph < 10.53:
                return ((6*data[2]["count"] + 10.53*data[3]["count"]) / (data[2]["count"] + data[3]["count"]))
            elif 10.53 < ph < 12.48:
                return ((10.53*data[3]["count"] + 12.48*data[4]["count"]) / (data[3]["count"] + data[4]["count"]))


print("Isoelectric point of the given protien is :", cal_iso_pt())
# The answers obtained is 4.65, it shows that the iso-electronic point of our protien.
# The ans may not be accurate as we have calculated the wieghted mean an did'nt knew the exact algorithm
