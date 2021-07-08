import os
import re
from json import dump, load
from bs4 import BeautifulSoup
from openbabel.pybel import readstring
from selfies import encoder

# from pymongo.errors import DuplicateKeyError, ConnectionFailure
# from atomate.qchem.database import QChemCalcDb


class ReaxysParser2:
    def __init__(self, base_dir):
        """
        ReaxysScraper
        :param base_dir: Base directory for parsing and generating data.
        """

        self.base_dir = base_dir

    def combiner(self, file_list, new_file_name):
        """
        Combines multiple JSON files into one file and saves to new JSON file
        :param file_list: JSON file list to be combined
        :param new_file_name: JSON file name for new combined file
        :return: dict of all combined JSON files
        """

        return_list = []

        for file in file_list:

            with open(file+".json") as f:

                rn = load(f)
                return_list.extend(rn)

        # saves to JSON and also returns the dict
        self.save_to_json(return_list, new_file_name)
        return return_list

    def replace_rn(self, json_file, substance_json):
        """
        Adds SMILES and SELFIES strings for reagents, catalysts, and solvents in JSON file
        :param json_file: JSON file produced from first parser script
        :param substance_json: JSON file of Bernstein Registry numbers and corresponding SMILES strings
        :return: list of reaction dicts including all SMILES and SELFIES
        """

        with open(substance_json) as f:
            rn = load(f)
        with open(json_file) as f:
            reactions = load(f)

        removed_reac_indices = []

        # iterates through every reaction
        # adds SMILES and SELFIES using Beilstein Registry numbers to every procedure
        # procedures with indices that do not correspond to SMILES or SELFIES are deleted
        # indices of reactions with no procedures are collected to be deleted
        for reaction_index, reaction in enumerate(reactions):

            removed_proc_indices = []

            # iterates through every procedure in each reaction
            # adds SMILES and SELFIES using Beilstein Registry numbers for catalysts, solvents, and reagents
            # collects procedures indices with lacking registry numbers to be deleted later
            for procedure_index, procedure in enumerate(reaction["procedures"]):

                # error count increases if a registry number does not correspond to a SMILES string
                error_cnt = 0

                reagents = procedure["reagents"]
                catalysts = procedure["catalysts"]
                solvents = procedure["solvents"]

                reagents_smi = []
                reagents_slf = []
                for reagent_rn in reagents["beilstein_rn"]:
                    try:
                        # tries to collect SMILES string and convert to SELFIES
                        reagent_smi = rn[str(reagent_rn)]
                        reagents_smi.append(reagent_smi)
                        reagents_slf.append(encoder(reagent_smi))
                    except:
                        error_cnt = error_cnt + 1
                # adds SELFIES and SMILES to reagents dictionary
                reaction["procedures"][procedure_index]["reagents"]["smiles"] = reagents_smi
                reaction["procedures"][procedure_index]["reagents"]["selfies"] = reagents_slf

                catalysts_smi = []
                catalysts_slf = []
                for catalyst_rn in catalysts["beilstein_rn"]:
                    try:
                        # tries to collect SMILES string and convert to SELFIES
                        catalyst_smi = rn[str(catalyst_rn)]
                        catalysts_smi.append(catalyst_smi)
                        catalysts_slf.append(encoder(catalyst_smi))
                    except:
                        error_cnt = error_cnt + 1
                # adds SELFIES and SMILES to catalysts dictionary
                reaction["procedures"][procedure_index]["catalysts"]["smiles"] = catalysts_smi
                reaction["procedures"][procedure_index]["catalysts"]["selfies"] = catalysts_slf

                solvents_smi = []
                solvents_slf = []
                for solvent_rn in solvents["beilstein_rn"]:
                    try:
                        # tries to collect SMILES string and convert to SELFIES
                        solvent_smi = rn[str(solvent_rn)]
                        solvents_smi.append(solvent_smi)
                        solvents_slf.append(encoder(solvent_smi))
                    except:
                        error_cnt = error_cnt + 1
                # adds SELFIES and SMILES to solvents dictionary
                reaction["procedures"][procedure_index]["solvents"]["smiles"] = solvents_smi
                reaction["procedures"][procedure_index]["solvents"]["selfies"] = solvents_slf

                if error_cnt != 0:
                    removed_proc_indices.append(procedure_index)

            # reverse sorts procedure indices and deletes from the back
            removed_proc_indices.sort(reverse=True)
            for removed_proc_index in removed_proc_indices:
                reaction["procedures"].pop(removed_proc_index)

            # sets changes back into the original reactions list
            reactions[reaction_index] = reaction

            # if zero procedures are left after this, adds index to list of reactions to be deleted
            if len(reaction["procedures"]) == 0:
                removed_reac_indices.append(reaction_index)

        # reverse sorts reaction indices and deletes from the back
        removed_reac_indices.sort(reverse=True)
        for removed_reac_index in removed_reac_indices:
            reactions.pop(removed_reac_index)

        return reactions

    def delete_duplicates(self, json_file):
        """
        Given a JSON file of new reactions, checks for duplicate reactions and procedures and returns a list of dicts
        :param json_file: JSON file of new reactions to be checked for duplicates
        :return: list of reaction dicts with no duplicate procedures or reactions
        """
        with open(json_file) as f:
            new_rxns = load(f)

        total_reactions = []

        # adds the first new reaction to total_reactions list
        if len(new_rxns) !=0:
            total_reactions = [new_rxns[0]]
            new_rxns.pop(0)

        # iterates through every new reaction
        # adds the procedures of the non-unique reaction to the first instance of the reaction
        # collects non-unique reaction indices to delete afterward
        for new_rxn in new_rxns:

            # a nonzero count will indicate that a reaction is a duplicate reaction
            count = 0

            # iterates through the old reactions to check if the new reaction matches any of the old reactions
            for old_rxn_index, old_rxn in enumerate(total_reactions):

                # checks if reaction ID is the same
                if old_rxn["reaction_id"] == new_rxn["reaction_id"]:
                    total_reactions[old_rxn_index]["procedures"].extend(new_rxn["procedures"])
                    count+=1
                    continue

                # checks if reactant ID is not the same but if reactants and products are the same
                elif set(old_rxn["reactants"]["beilstein_rn"]) == set(new_rxn["reactants"]["beilstein_rn"]) and set(old_rxn["products"]["beilstein_rn"]) == set(new_rxn["products"]["beilstein_rn"]):
                    total_reactions[old_rxn_index]["procedures"].extend(new_rxn["procedures"])
                    rxn_id = [old_rxn["reaction_id"],new_rxn["reaction_id"]]
                    total_reactions[old_rxn_index]["reaction_id"] = rxn_id
                    count+=1
                    continue

            # if count does not increase (no duplicates), new reaction is added to list of total reactions
            if count == 0:
                total_reactions.append(new_rxn)

        # iterates through total reactions to look for duplicate procedures within reactions
        for rxn in total_reactions:

            # iterates through all procedures
            for procedure_index, procedure in enumerate(rxn["procedures"]):

                # for each procedure, checks through all other procedures if any are duplicates
                # collects the index of any duplicates to delete afterward
                removed_proc_indices = []
                for test_procedure_index, test_procedure in enumerate(rxn["procedures"]):

                    # ensures the procedures are not equal to each other
                    if procedure_index != test_procedure_index:

                        # if reagents, solvents, and experiment text are equal, collects index of test procedure
                        if set(procedure["reagents"]["beilstein_rn"]) == set(test_procedure["reagents"]["beilstein_rn"]) and set(procedure["solvents"]["beilstein_rn"]) == set(test_procedure["solvents"]["beilstein_rn"]) and procedure["exp_text"] == test_procedure["exp_text"]:
                            removed_proc_indices.append(test_procedure_index)

                # reverse sorts procedure indices and deletes from the back
                removed_proc_indices.sort(reverse=True)
                for removed_proc_index in removed_proc_indices:
                    rxn["procedures"].pop(removed_proc_index)

        return total_reactions


    def get_balanced_reactions(self, reaction_json):
        with open(reaction_json) as f:
            reactions = load(f)

        return reactions


    def save_to_json(self, dict_list, filename):
        """
        Saves list of dicts to a JSON file
        :param dict_list: list of dicts to be written to JSON file.
        :param filename: name of JSON file without extension
        :return: bool if successful in saving file
        """
        try:
            with open(filename + ".json", "w") as outfile:
                dump(dict_list, outfile)
            return True
        except:
            return False


# start of code

parser = ReaxysParser2("/Users/nikitaredkar/Documents/COLLEGE/Research/PerssonLab/DATA/")
filenames = ["oxidation_1_2640", "reduction_1_2524", "electrolysis_1_4911", "electrochemical_85001_85072"]

i = 0
while i < 80000:
    filename = "electrochemical_" + str(i+1) + "_" + str(i+5000)
    filenames.append(filename)
    i+=5000

for index, file in enumerate(filenames):
    results = parser.replace_rn(file+".json", "reagents.json")
    parser.save_to_json(results, file+"_updated")
    filenames[index] = file+"_updated"

parser.combiner(filenames, "total_reactions_with_duplicates")

results = parser.delete_duplicates("total_reactions_with_duplicates.json")
parser.save_to_json(results, "total_reactions_without_duplicates")









