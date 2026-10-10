"""Reproduce audit and exploratory runs; no email, upload or publication."""
import subprocess
import sys
from pipeline import ROOT


def main():
    commands=[['acquire.py'],['fetch_damage.py'],['pipeline.py','inventory'],['audit_damage.py']]
    for h in [3,4]:
        commands.append(['pipeline.py','analyse','data/private/buildings.csv','data/private/hypothesis_assignments.csv',
                         f'data/private/damage_map_heavy{h}.csv',f'data/private/exploratory_heavy{h}'])
    commands.append(['summarize.py'])
    for args in commands:subprocess.run([sys.executable,*args],cwd=ROOT,check=True)


if __name__=='__main__':main()
