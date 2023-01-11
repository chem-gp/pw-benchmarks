from tblite.ase import TBLite
from ase.atoms import Atoms
from ase.io import write, read
import numpy as np
import ase.units

# atoms = Atoms(
#     symbols="C4O8",
#     positions=np.array(
#         [
#             [0.9441259872, 0.9437851680, 0.9543505632],
#             [3.7179966528, 0.9556570368, 3.7316862240],
#             [3.7159517376, 3.7149292800, 0.9692330016],
#             [0.9529872864, 3.7220864832, 3.7296981120],
#             [1.6213905408, 1.6190616096, 1.6313879040],
#             [0.2656685664, 0.2694175776, 0.2776540416],
#             [4.3914553920, 1.6346256864, 3.0545920000],
#             [3.0440834880, 0.2764611744, 4.4080419264],
#             [4.3910577696, 3.0416409504, 0.2881058304],
#             [3.0399936576, 4.3879335936, 1.6497353376],
#             [0.2741322432, 4.4003734944, 3.0573754368],
#             [1.6312174944, 3.0434586528, 4.4023048032],
#         ]
#     ),
#     cell=np.array([5.68032, 5.68032, 5.68032]),
#     pbc=np.array([True, True, True]),
# )

# crystal=read("../structures/HIK-143 293K-activated.cif")
crystal=read("../structures/KTU-183_2_auto.cif")


# crystal=read("../structures/HIK-143 293K-activated.cif")
atoms = crystal
atoms.calc = TBLite(method="GFN2-xTB", verbosity=0)
atoms.set_pbc(True)

print( atoms.get_forces()  )