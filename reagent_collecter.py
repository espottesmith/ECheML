import os
import re
from json import dump, load
from bs4 import BeautifulSoup
from openbabel.pybel import readstring
from selfies import encoder

class ReagentCollecter:
    def __init__(self, base_dir):
        """
        ReagentCollector
        :param base_dir: Base directory for parsing and generating data.
        """

        self.base_dir = base_dir

    def reagent_parser(self, filename, initial_dictlist = dict()):
        """
        Parses Beilstein Registry Number data to output a list of dictionaries relating number and SMILES
        :param filename: name of XML file with extension
        :param dict_list: list of dicts to be written to JSON file.
        :return: dict with Registry Number/SMILES string pairing
        """
        results = initial_dictlist

        filepath = os.path.join(self.base_dir, filename)

        with open(filepath, 'r') as fileobj:
            xml = fileobj.read()
            parsed = BeautifulSoup(xml, "lxml-xml")

            substances = parsed.find_all("substance")

            # iterates for every substance in file
            for substance in substances:

                # pulls the Beilstein Registry number
                rn = substance.find("IDE.XRN").text

                # finds all molfiles, indicated by "YY.STR"
                molfiles = substance.find_all("YY.STR")

                # V2000 files always follow V3000 files if V3000 files are available
                # if V3000 files are not available, V2000 files are the first and only file
                # we are only using V2000 files
                try:

                    # if possible, converts the second file (the V2000 file) to SMILES
                    mymol = readstring("mol", molfiles[1].text)
                    smi = mymol.write(str("smi")).split()[0]

                    results[rn] = smi
                except:

                    # occurs when only zero or one file listed for specific registry number
                    try:

                        # if possible, converts the first file (the V2000 file) to SMILES
                        mymol = readstring("mol", molfiles[0].text)
                        smi = mymol.write(str("smi")).split()[0]

                        results[rn] = smi
                    except:

                        # if not file is listed for the registry number, skips that number
                        continue

        return results

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


# start of parsing code

parser = ReagentCollecter("/Users/nikitaredkar/Documents/COLLEGE/Research/PerssonLab/DATA/")
filenames = ["reagents_1_5000", "reagents_5001_10000", "reagents_10001_15000", "reagents_15001_end"]
results = dict()

for file in filenames:
    results = parser.reagent_parser(file+".xml", results)

parser.save_to_json(results, "reagents")
