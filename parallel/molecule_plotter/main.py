# Author: Alex
import sys
import os
from ase import io

if len(sys.argv) == 1:
    filename = "trajectory.gif"
elif sys.argv[1][-4:] != ".gif":
    filename = f"{sys.argv[1]}.gif"
else:
    filename = sys.argv[1]

if not os.path.exists("./results/plots"):
    os.makedirs("./results/plots")

try:
    atoms = io.read("./results/trajectory.xyz", index=":")
    io.write(f"./results/plots/{filename}", atoms, interval=50)
    print(f"{filename} was created")

except:
    print(
        "Structure trajectory could not be found. Possible reasons:\n - Did you try to plot it before creating it?\n - Did you move or rename the results/trajectory.xyz file?"
    )
