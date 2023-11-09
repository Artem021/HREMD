import os, subprocess, re

LMP = '/home/artem/LAMMPS/lammps-static/bin/lmp'

os.chdir('/home/artem/LAMMPS/test/PAM')

CMD_LINE = [LMP,'-i', 'PAM.lmp']


popen = subprocess.Popen(CMD_LINE, stdout=subprocess.PIPE, universal_newlines=True,stderr=None)
pat = re.compile('Time\s*Temp\s*PotEng\s*KinEng\s*Press\s*Volume\s*Density')
header_found = False
dE = []
dt = []
for stdout_line in iter(popen.stdout.readline, ""):
    if re.search(pat,stdout_line) != None:
        header_found = True
        continue
    if re.search('Performance',stdout_line) !=None and header_found:
        break
    if header_found:
        try:
            t, _, E, _, _, _, _ = [int(stdout_line.split()[0])] + [float(i) for i in stdout_line.split()[1:]]
            dE.append(E)
            dt.append(t)
        except:
            continue
        print(stdout_line)

return_code = popen.wait()
popen.stdout.close()



def write_input_lammps():
    pass