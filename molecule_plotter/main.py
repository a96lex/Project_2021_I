import os
from ase import io
import imageio

atoms = io.read('../results/trajectory.xyz', index=":")

images=[]
for i in range(1000):
    filename="p{}.png".format(i)
    io.write(filename,atoms[i])
    images.append(imageio.imread(filename))
    os.remove(filename)

imageio.mimsave('movie.gif', images)


