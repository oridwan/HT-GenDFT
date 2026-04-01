#!/usr/bin/env python
import os
import time
import numpy as np
import codecs
from pymatgen.core.periodic_table import Element

def wait_for_file(filepath, timeout=30, check_interval=0.2, expected_lines=None):
    """
    Wait for a file to exist and be readable on network filesystems.
    
    Increased defaults for HPC environments with heavy parallel workloads.
    
    Args:
        filepath: Path to file to wait for
        timeout: Maximum time to wait in seconds (default: 30s for HPC)
        check_interval: Time between checks in seconds (default: 0.2s)
        expected_lines: If provided, validate file has this many lines (for text files)
    
    Returns:
        bool: True if file exists and is valid, False if timeout
    """
    elapsed = 0
    while elapsed < timeout:
        if os.path.exists(filepath):
            try:
                # Use os.stat to force metadata cache refresh
                stat_info = os.stat(filepath)
                size = stat_info.st_size
                
                # Extra safety: verify file is not empty
                if size > 0:
                    # For PARCHG-m files with expected line count, validate completeness
                    if expected_lines is not None and filepath.endswith('-m'):
                        try:
                            with open(filepath, 'r') as f:
                                actual_lines = sum(1 for _ in f)
                            if actual_lines == expected_lines:
                                return True
                            # File exists but incomplete, keep waiting
                        except:
                            pass
                    else:
                        return True
            except:
                pass
        time.sleep(check_interval)
        elapsed += check_interval
    return False


class electride:
    """
    Parse electride from a directory containing ELFCAR, PARCHG and vasprun.xml
    """

    def __init__(self, path = './',
        vaspxml = 'vasprun.xml',
        cmd = 'bader ',
        ELF_min = 0.60):

        self.label = ['e0025', 'e05', 'e10', 'band0', 'band1']
        self.volume = [0]*len(self.label)
        self.error = None
        self.ELF_maxima = []

        # Save current directory and change to path for all file operations
        original_dir = os.getcwd()
        path = os.path.abspath(path)
        
        if not os.path.exists(os.path.join(path, 'ELFCAR')):
            self.error = 'ELFCAR not exist'
            return
        
        # Change to structure directory for all file operations
        os.chdir(path)
        
        try:
            ELFCAR_path = 'ELFCAR'
            PARCHG_path = 'PARCHG'
            clean_cmd1 = 'rm -f ELFCAR-m bader_log ACF.dat BCF.dat AVF.dat'
            clean_cmd2 = 'rm -f PARCHG-m bader_log ACF.dat BCF.dat AVF.dat'

            cell, coor, radii, grid = self.Read_ELFCAR(ELFCAR_path, ELF_min)
            
            # Wait for ELFCAR-m to be visible on network filesystem
            if not wait_for_file('ELFCAR-m', timeout=30):
                self.error = 'ELFCAR-m not created'
                return
            
            # Additional sync delay before running bader under heavy parallel load
            time.sleep(0.5)
            ret = os.system(cmd + 'ELFCAR-m > bader_log 2>&1')
            
            # Wait for BCF.dat with retry for network filesystem sync
            timeout = 30 if ret == 0 else 5
            if not wait_for_file('BCF.dat', timeout=timeout):
                if ret != 0:
                    self.error = f'bader command failed (exit code {ret >> 8}) for ELFCAR'
                else:
                    self.error = 'bader error in parsing ELFCAR'
            else:
                self.ELF_maxima, self.fac1 = self.Read_BCF('BCF.dat', radii, cell)
            os.system(clean_cmd1)

            if len(self.ELF_maxima) > 0:
                self.cell = cell
                for i, label in enumerate(self.label):
                    fpath = PARCHG_path + '-' + label
                    if os.path.exists(fpath):
                        expected_lines = self.ModifyPARCHG(self.ELF_maxima, fpath)
                        
                        # Wait for PARCHG-m to be fully written and validated
                        if not wait_for_file('PARCHG-m', timeout=30, expected_lines=expected_lines):
                            self.error = f'PARCHG-m not created or incomplete for {fpath}'
                            continue
                        
                        # Additional sync delay before running bader under heavy parallel load
                        time.sleep(0.5)
                        ret = os.system(cmd + 'PARCHG-m -vac off > bader_log 2>&1')

                        # Wait for ACF.dat with retry for network filesystem sync
                        # If bader failed (non-zero return), don't wait as long
                        timeout = 30 if ret == 0 else 5
                        if wait_for_file('ACF.dat', timeout=timeout):
                            charge = self.Read_ACF('ACF.dat')
                            os.system(clean_cmd2)
                            self.volume[i] = self.parse_charge(charge)
                        else:
                            if ret != 0:
                                self.error = f'Bader command failed (exit code {ret >> 8}) for {fpath}'
                            else:
                                self.error = f'Bader analysis failed for {fpath}'
        finally:
            # Always return to original directory
            os.chdir(original_dir)

    def parse_charge(self, charge):
        totalv = 0.0
        for i in range(len(self.ELF_maxima)):
            if charge[i][0] > 0.01:
                totalv += charge[i][2]
        volume = totalv/np.linalg.det(self.cell)*100
        return volume

    @staticmethod
    def Read_ACF(filename):
        with open(filename, 'rb') as f:
            input_content = f.readlines()
        chg = []
        for i in range(len(input_content)):
            s = input_content[i].split()
            if s[0].isdigit():
               a=[float(f) for f in s]
               chg.append(a[4:])
        return np.array(chg)

    @staticmethod
    def Read_BCF(filename, radii, cell):
        with open(filename, 'rb') as f:
            input_content = f.readlines()
        pos = []
        count = 0
        cell = np.linalg.inv(cell)
        for i in range(len(input_content)):
            s = input_content[i].split()
            if s[0].isdigit():
               a=[float(f) for f in s]
               if 1.2*radii[int(a[-2])-1] < a[-1]:
                  count += 1
                  if len(pos) < 4:
                     pos.append(np.dot(np.array(a[1:4]), cell))
        if len(pos) > 0:
           fac = count/len(pos)
        else:
           fac = 1.0
        return np.array(pos), fac

    @staticmethod
    def ModifyPARCHG(pos, filename):
        """
        This module modifies PARCHG for bader analysis.
        Returns the expected line count for validation.
        """
        with codecs.open(filename, 'r', encoding='utf-8') as f1:
            input_content = f1.readlines()
        
        expected_lines = len(input_content) + len(pos)
        
        with open('PARCHG-m', 'w') as f2:
            f2.writelines(input_content[:5])
            f2.write(' H '+input_content[5])
            f2.write(str(len(pos)) + ' ' + input_content[6])
            f2.write(input_content[7])
            for coor in pos:
                f2.write('%10.6f %9.6f %9.6f\n' % (coor[0], coor[1], coor[2]))
            f2.writelines(input_content[8:])
            f2.flush()
            os.fsync(f2.fileno())
        
        return expected_lines

    @staticmethod
    def Read_ELFCAR(filename, ELF_min):
        """
        This module reads ELFCAR
        """
        with open(filename, 'rb') as f:
            input_content = f.readlines()
        
        count = 0
        cell = []
        coor = []
        ELF_raw = []
        N_atoms = 0
        grid = []
        
        with open('ELFCAR-m', 'w') as f1:
            for line in input_content:
                line=str(line,'utf-8')
                count = count + 1
                if count < 3:
                   f1.write(line)
                elif (count>=3) and (count<=5):
                   cell.append([float(f) for f in line.split()])
                   f1.write(line)
                elif count==6:
                   f1.write(line)
                   symbol = line.split()
                elif count==7:
                   f1.write(line)
                   numIons = [int(f) for f in line.split()]
                   N_atoms = sum(numIons)
                   f1.write('Direct\n')
                elif (count>=9) and (count<9+N_atoms):
                   f1.write(line)
                   coor.append([float(f) for f in line.split()])
                elif count == 10+N_atoms:
                   f1.write('\n')
                   f1.write(line)
                   grid = [int(f) for f in line.split()]
                elif count > 10+N_atoms:
                   ELF_raw = line.split()
                   for i, f in enumerate(ELF_raw):
                       if float(f)<ELF_min:
                          f = '0.00000'
                       if i==0:
                          f1.write('%8s' % (f))
                       else:
                          f1.write('%12s' % (f))
                   f1.write('\n')
            f1.flush()
            os.fsync(f1.fileno())

        radii = np.array([])
        for ele in range(len(symbol)):
            # Use atomic_radius from pymatgen (in Angstroms)
            radius = Element(symbol[ele]).atomic_radius
            if radius is None:
                radius = 1.0  # Default fallback for elements without atomic_radius
            else:
                radius = float(radius)  # Convert FloatWithUnit to plain float
            tmp = radius * np.ones(numIons[ele])
            radii = np.append(radii, tmp)
        cell = np.array(cell)
        return cell, coor, radii, grid

def get_subdir(a_dir):
   return sorted([name for name in os.listdir(a_dir)
           if os.path.isdir(os.path.join(a_dir, name))
           and name not in ['mp_mattersim_cache', '.git', '__pycache__']])

if __name__ == "__main__":
    from optparse import OptionParser
    import pandas as pd
    from tabulate import tabulate

    parser = OptionParser()
    parser.add_option("-d", "--directory", dest="dir", default='./',
                      help="directory containing PARCHG and ELFCAR", metavar="dir")
    parser.add_option("-p",  "--pdir", dest='pdir',
                  help="by parent directory")
    parser.add_option("-b", "--bader-exe", dest='bader_exe', default='bader',
                  help="path to bader executable (default: bader)")
    parser.add_option("-t", "--threshold", dest='threshold', type='float', default=0.6,
                  help="ELF threshold value (default: 0.6)")

    (options, args) = parser.parse_args()
    total_dir = []
    if options.pdir is None:
       total_dir.append(options.dir)
    else:
       total_dir =  get_subdir(options.pdir)
       os.chdir(options.pdir)

    col_name = {'formula':[]}
    for label in ['e0025','e05', 'e10', 'band0', 'band1']:
        col_name[label] = []

    for subdir in total_dir:
        print(subdir)
        test = electride(subdir, cmd=options.bader_exe + ' ', ELF_min=options.threshold)
        
        # Extract structure name from path (e.g., "Li10B1N4_s001" from full path)
        # Handles both "path/to/Li10B1N4_s001/ELF" and "Li10B1N4_s001"
        import re
        struct_name = subdir
        match = re.search(r'([A-Z][a-z]?\d+(?:[A-Z][a-z]?\d+)*_s\d+)', subdir)
        if match:
            struct_name = match.group(1)
        else:
            # Fallback: just use basename
            struct_name = os.path.basename(subdir.rstrip('/'))
        
        col_name['formula'].append(struct_name)
        
        # Print error if any
        if test.error:
            print(f"  Warning: {test.error}")

        for i, label in enumerate(['e0025','e05', 'e10', 'band0', 'band1']):
            col_name[label].append(test.volume[i])

    df = pd.DataFrame(col_name)
    print(tabulate(df, headers='keys', tablefmt='psql'))
    df.to_csv('electride_analysis.csv')
