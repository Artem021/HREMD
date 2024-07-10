from engines import runMD
import re,os,shutil

arg = ('/home/artem/LAMMPS_TEST/macro/10.07/0.000e+00', 'lammps.inp', '/home/artem/LAMMPS_TEST/macro.data', 'lammps.error', 'traj.xyz', 10, re.compile('Time\\s+Temp\\s+PotEng\\s+KinEng\\s+Press\\s+Volume\\s+Density'), 2)

res = runMD(arg)