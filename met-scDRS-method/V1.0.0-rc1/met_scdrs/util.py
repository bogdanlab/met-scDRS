def convert_species_name(species):
    if species in ["Mouse", "mouse", "Mus_musculus", "mus_musculus", "mmusculus"]:
        return "mmusculus"
    if species in ["Human", "human", "Homo_sapiens", "homo_sapiens", "hsapiens"]:
        return "hsapiens"
    raise ValueError("species name '%s' is not supported" % species)