# Author: Alex
import sys
import os
from ase import io

if len(sys.argv) == 1:
    filename = "trajectory.gif"
elif sys.argv[1][-4:] != ".gif":
    filename = "{}.gif".format(sys.argv[1])
else:
    filename = sys.argv[1]

if not os.path.exists("./results/plots"):
    os.makedirs("./results/plots")

try:
    atoms = io.read("./results/trajectory.xyz", index=":")
    io.write("./results/plots/{}".format(filename), atoms, interval=50)
    print("{} was created".format(filename))
except:
    print(
        "Structure trajectory could not be found. Possible reasons:\n - Did you try to plot it before creating it?\n - Did you move or rename the results/trajectory.xyz file?"
    )
