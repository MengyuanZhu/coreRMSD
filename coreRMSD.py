import openbabel
from rdkit import Chem
from rdkit.Chem import AllChem
import math


def getRMSD(distances):  #get final rmsd value
	numDistances=len(distances)
	itrDistance=0
	deviations=0
	while itrDistance<numDistances:
		deviations=deviations+distances[itrDistance]
		itrDistance=itrDistance+1
	return round(math.sqrt(deviations/numDistances),3)

def getDistance(pos1,pos2):  #calculate the distance square of two atoms
	return (pos1.x-pos2.x)**2+(pos1.y-pos2.y)**2+(pos1.z-pos2.z)**2

obConversion=openbabel.OBConversion()
obConversion.SetInFormat("mol2")
print sb
obmol=openbabel.OBMol()
notatend=obConversion.ReadFile(obmol,"csr_ladi-results.mol2")
while notatend:
	#rdmol=openbabel.OBMolToRWMol(obmol)

	obmol=openbabel.OBMol()
	notatend=obConversion.Read(obmol)


m1=Chem.MolFromMol2File("test1.mol2")
m2=Chem.MolFromMol2File("test2.mol2")
core=Chem.MolFromMol2File("core.mol2")
core1=Chem.ReplaceSidechains(m1,core)
core2=Chem.ReplaceSidechains(m2,core)
atoms1=core1.GetAtoms()
atoms2=core2.GetAtoms()

i=0
myarray=[]
while i < len(atoms1):
	pos1=core1.GetConformer().GetAtomPosition(i)
	pos2=core2.GetConformer().GetAtomPosition(i)
	myarray.append(getDistance(pos1,pos2))	
	i=i+1
print getRMSD(myarray)






