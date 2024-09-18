# mutate.py

from pymol import cmd

def mutate(molecule, chain, resi, target="CYS", mutframe="1"):
    target = target.upper()
    cmd.fetch(molecule)
    cmd.wizard("mutagenesis")
    cmd.do("refresh_wizard")
    cmd.get_wizard().set_mode(target)
    selection = f"/{molecule}//{chain}/{resi}"
    cmd.get_wizard().do_select(selection)
    cmd.frame(str(mutframe))
    cmd.get_wizard().apply()
    cmd.set_wizard()
    cmd.save(f"{molecule}_{resi}_{target}.pdb", molecule)
cmd.extend("mutate", mutate)

