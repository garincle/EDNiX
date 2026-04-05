#import
import os
import xml.etree.ElementTree as ET

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile


def read(node, indent=0):
    print("  " * indent + f"{node.tag} {node.attrib} {node.text}")
    for child in node:
        read(child, indent + 1)

def renamefiles(specfile,oldname,newname):

    tree = ET.parse(specfile)
    root = tree.getroot()

    for datafile in root.findall(".//DataFile"):
        if datafile.text and oldname in datafile.text:
            datafile.text = datafile.text.replace(oldname, newname)

    tree.write(specfile, encoding="UTF-8", xml_declaration=True)

def create(listfile,specfile):

    # make sure all the elements of the list can be found with a relative path
    # from the localisation of the specfile

    # Create root element
    root = ET.Element("CaretSpecFile", Version="1.0")

    # Add MetaData section
    metadata = ET.SubElement(root, "MetaData")

    for ftype in listfile:

        fname = ftype.split(".")

        if fname[1].low() == 'l':
            structname = 'CortexLeft'
        elif fname[1].low() == 'r':
            structname = 'CortexLeft'
        else:
            structname = 'Invalid'

        if fname[-1] == 'nii':
            if fname[-2] == 'dlabel':
                structname = 'All'
                Dataname = 'CONNECTIVITY_DENSE_LABEL'
            elif fname[-2] == 'dscalar':
                structname = 'All'
                Dataname = 'CONNECTIVITY_DENSE_SCALAR'

        elif fname[-1] == 'gii':
            if fname[-2] == 'surf':
                Dataname = 'SURFACE'
            elif fname[-2] == 'shape' or fname[-2] == 'func':
                Dataname = 'METRIC'
            elif fname[-2] == 'label':
                Dataname = 'LABEL'

        elif fname[-1] == 'gz' and fname[-2] == 'nii':
            Dataname = 'VOLUME'

        datafile = ET.SubElement(
            root,
            "DataFile",
            Structure=structname,
            DataFileType=Dataname,
            Selected="true"
        )
        datafile.text = f"\n      {fname}\n   "

    # Create tree and write file
    tree = ET.ElementTree(root)
    tree.write(specfile, encoding="UTF-8", xml_declaration=True)

    print("Spec file created.")