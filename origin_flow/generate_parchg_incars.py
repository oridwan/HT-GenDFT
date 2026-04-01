#!/usr/bin/env python3
"""
Generate PARCHG INCAR files from SC vasprun.xml

This helper script is called by the SPE SLURM job after SC calculation completes.
It parses vasprun.xml-SC to extract NELECT and grid parameters, then generates
all 5 PARCHG INCAR files (band0, band1, e0025, e05, e10).

Usage:
    python3 generate_parchg_incars.py <vasprun_xml_path> <reciprocal_density>

Example:
    python3 generate_parchg_incars.py vasprun.xml-SC 64
    python3 generate_parchg_incars.py vasprun.xml-SC 250  # For refined workflow
"""

import sys
from pathlib import Path
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.sets import MPStaticSet


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 generate_parchg_incars.py <vasprun_xml_path> <reciprocal_density>")
        sys.exit(1)
    
    vasprun_path = sys.argv[1]
    reciprocal_density = int(sys.argv[2])
    
    if not Path(vasprun_path).exists():
        print(f"ERROR: {vasprun_path} not found")
        sys.exit(1)
    
    try:
        # Parse vasprun.xml-SC to get parameters
        vr = Vasprun(vasprun_path, parse_dos=False, parse_eigen=True)
        structure = vr.final_structure
        
        nelect = int(vr.parameters['NELECT'])
        is_soc = vr.parameters.get('LSORBIT', False)
        fac = 1 if is_soc else 2
        
        if nelect % 2 == 0:
            iband = int(nelect / fac)
        else:
            iband = int(nelect / fac) + 1
        
        grid = [vr.parameters['NGXF'], vr.parameters['NGYF'], vr.parameters['NGZF']]
        
        print(f"NELECT={nelect}, VBM band index={iband}, Grid={grid}")
        
        # Keep PARCHG-specific controls, but use the same smearing policy as the SC step.
        base_settings = {
            'PREC': 'Accurate',
            'ALGO': 'Normal',
            'EDIFF': 1e-6,
            'IBRION': -1,
            'NSW': 0,
            'ISMEAR': 0,
            'SIGMA': 0.05,
            'ISPIN': 1,
            'ISYM': 0,
            'NPAR': 4,
            'LAECHG': False,
            'LASPH': True,
            'LVHAR': False,
            'LCHARG': False,
            'LWAVE': False,
            'ISTART': 1,
            'ICHARG': 11,
        }
        
        # Generate PARCHG INCARs for band-specific
        for label, iband_offset in [('band0', 0), ('band1', -1)]:
            settings = base_settings.copy()
            settings.update({
                'IBAND': iband + iband_offset,
                'LPARD': True,
                'NGXF': grid[0],
                'NGYF': grid[1],
                'NGZF': grid[2],
            })
            vis = MPStaticSet(structure, user_incar_settings=settings, 
                             user_kpoints_settings={'reciprocal_density': reciprocal_density})
            vis.incar.write_file(f'INCAR-PARCHG-{label}')
            print(f"Generated INCAR-PARCHG-{label}")
        
        # Generate PARCHG INCARs for energy windows
        for label, e_val in [('e0025', 0.025), ('e05', 0.5), ('e10', 1.0)]:
            settings = base_settings.copy()
            settings.update({
                'EINT': f'{-e_val} 0.025',
                'NBMOD': -3,
                'LPARD': True,
                'NGXF': grid[0],
                'NGYF': grid[1],
                'NGZF': grid[2],
            })
            vis = MPStaticSet(structure, user_incar_settings=settings,
                             user_kpoints_settings={'reciprocal_density': reciprocal_density})
            vis.incar.write_file(f'INCAR-PARCHG-{label}')
            print(f"Generated INCAR-PARCHG-{label}")
        
        print("All PARCHG INCARs generated successfully")
        
    except Exception as e:
        print(f"ERROR generating PARCHG INCARs: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
