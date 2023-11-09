import pyxyz
import pandas as pd
fxyz = "/home/artem/HREMD_algorithm/cp2k/rmsd/opt/vel_soft-8-xtbopt/structures-opt-gfnff-filter.xyz"

p = pyxyz.Confpool()
p.include_from_file(fxyz)
print(f"Loaded {p.size} conformers")

p['A'] = lambda m: m.v(53,58,112)
p['D'] = lambda m: m.z(50,14,15,16)

p.sort('A')

p.descr = lambda m: f'A = {m["A"]}; D = {m["D"]}'

df = p.as_table()
df = pd.DataFrame(df)
print(df)

# p.save('/home/artem/HREMD_algorithm/cp2k/rmsd/crest_results/structures-opt-gfnff-filter-Vang-Z.xyz')