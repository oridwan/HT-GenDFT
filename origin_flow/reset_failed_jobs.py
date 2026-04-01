#!/usr/bin/env python3
"""
Reset failed VASP jobs to allow retry.

This script resets structures that failed at various stages:
- RELAX_FAILED -> PENDING (retry from scratch)
- SC_FAILED -> RELAX_DONE (retry SC stage in SPE)
- PARCHG_FAILED -> SC_DONE (retry PARCHG stage in SPE)
- ELF_FAILED -> PARCHG_DONE (retry ELF stage in SPE)

Usage:
    python3 reset_failed_jobs.py [--db workflow.json] [--stage STAGE] [--clean] [--dry-run]

Options:
    --db PATH           Path to workflow.json (default: ./VASP_JOBS/workflow.json)
    --stage STAGE       Only reset specific stage: RELAX, SC, PARCHG, or ELF (default: all)
    --clean             Clean calculation directories to proper restart state:
                        - RELAX_FAILED: Remove entire Relax directory
                        - SC_FAILED: Remove all SPE directory contents (keeps Relax/CONTCAR)
                        - PARCHG_FAILED: Remove PARCHG and intermediate files (keeps SC outputs)
                        - ELF_FAILED: Remove ELF and intermediate files (keeps PARCHG files)
    --dry-run           Show what would be done without making changes
    --list              List all failed structures and exit

Examples:
    # Dry run to see what would be reset
    python3 reset_failed_jobs.py --dry-run

    # List all failed structures
    python3 reset_failed_jobs.py --list

    # Reset only RELAX_FAILED structures (no directory cleanup)
    python3 reset_failed_jobs.py --stage RELAX

    # Reset all failed jobs and clean directories to proper restart state
    python3 reset_failed_jobs.py --clean

    # Reset ELF failures with cleanup (removes old ELFCAR, keeps PARCHG-* files)
    python3 reset_failed_jobs.py --stage ELF --clean
"""

import json
import shutil
import argparse
from pathlib import Path
from collections import defaultdict


def list_failed_structures(data):
    """List all failed structures grouped by failure stage."""
    failed_states = ['RELAX_FAILED', 'SC_FAILED', 'PARCHG_FAILED', 'ELF_FAILED']
    
    failed_by_stage = defaultdict(list)
    for struct_id, sdata in data['structures'].items():
        if sdata['state'] in failed_states:
            failed_by_stage[sdata['state']].append({
                'id': struct_id,
                'composition': sdata.get('composition', 'N/A'),
                'error': sdata.get('error', 'No error message'),
                'dir': sdata.get(f"{sdata['state'].split('_')[0].lower()}_dir", 'N/A')
            })
    
    if not any(failed_by_stage.values()):
        print("No failed structures found!")
        return
    
    print("\n" + "="*80)
    print("Failed Structures Summary")
    print("="*80)
    
    for state in failed_states:
        if state in failed_by_stage:
            structures = failed_by_stage[state]
            print(f"\n{state}: {len(structures)} structures")
            print("-" * 80)
            for s in structures[:5]:  # Show first 5
                print(f"  {s['id']:<25} {s['composition']:<20}")
                if s['error'] and s['error'] != 'No error message':
                    error_short = s['error'][:60] + '...' if len(s['error']) > 60 else s['error']
                    print(f"    Error: {error_short}")
            if len(structures) > 5:
                print(f"  ... and {len(structures) - 5} more")
    
    print("\n" + "="*80)
    print(f"Total failed: {sum(len(v) for v in failed_by_stage.values())}")
    print("="*80 + "\n")


def reset_failed_jobs(db_path, stage_filter=None, clean=False, dry_run=False):
    """
    Reset failed jobs to allow retry.
    
    Args:
        db_path: Path to workflow.json
        stage_filter: Only reset specific stage ('RELAX', 'SC', 'PARCHG', or None for all)
        clean: Remove failed directories and error markers
        dry_run: Don't make actual changes
    """
    db_path = Path(db_path)
    
    if not db_path.exists():
        print(f"Error: Database not found at {db_path}")
        return
    
    print(f"Loading database: {db_path}")
    with open(db_path, 'r') as f:
        data = json.load(f)
    
    # Define reset mappings: failed_state -> (new_state, job_id_to_clear, dir_to_check)
    # Note: SC/PARCHG/ELF all reset to RELAX_DONE because SPE is a unified job
    # that runs all three stages sequentially. The workflow only triggers submit_spe()
    # when state is RELAX_DONE or RELAX_TMOUT.
    reset_map = {
        'RELAX_FAILED': ('PENDING', 'relax_job_id', 'relax_dir'),
        'SC_FAILED': ('RELAX_DONE', 'spe_job_id', 'spe_dir'),
        'PARCHG_FAILED': ('RELAX_DONE', 'spe_job_id', 'spe_dir'),
        'ELF_FAILED': ('RELAX_DONE', 'spe_job_id', 'spe_dir'),
    }
    
    # Filter by stage if requested
    if stage_filter:
        stage_upper = stage_filter.upper()
        reset_map = {k: v for k, v in reset_map.items() if k.startswith(stage_upper)}
        if not reset_map:
            print(f"Error: Invalid stage '{stage_filter}'. Must be RELAX, SC, PARCHG, or ELF")
            return
    
    reset_counts = defaultdict(int)
    cleaned_dirs = []
    
    print("\nScanning for failed structures...")
    if dry_run:
        print("DRY RUN MODE - no changes will be made\n")
    
    if not clean:
        print("NOTE: Running without --clean flag")
        print("      Job states will be reset but directories won't be cleaned")
        print("      Use --clean to remove old calculation files for fresh restart\n")
    
    for struct_id, sdata in data['structures'].items():
        state = sdata['state']
        
        if state not in reset_map:
            continue
        
        new_state, job_field, dir_field = reset_map[state]
        job_dir = Path(sdata[dir_field])
        
        print(f"\n{struct_id} ({state}):")
        print(f"  Current state: {state}")
        print(f"  Will reset to: {new_state}")
        print(f"  Directory: {job_dir}")
        
        if not dry_run:
            # Update state
            sdata['state'] = new_state
            
            # Clear job IDs
            sdata[job_field] = None
            
            # Clear error message
            sdata['error'] = None
            
            reset_counts[state] += 1
        else:
            reset_counts[state] += 1
        
        # Clean up directories if requested
        if clean and job_dir.exists():
            print(f"  Cleaning directory to restart from {new_state}...")
            
            if state == 'RELAX_FAILED':
                # Remove entire Relax directory to start fresh
                print(f"    - Removing entire Relax directory")
                if not dry_run:
                    shutil.rmtree(job_dir)
                    cleaned_dirs.append(str(job_dir))
            
            elif state == 'SC_FAILED':
                # Relax/CONTCAR is in separate relax_dir and will be preserved
                print(f"    - Removing all SPE directory contents for SC retry")
                relax_dir = Path(sdata['relax_dir'])
                print(f"    - Relax/CONTCAR will be preserved at: {relax_dir / 'CONTCAR'}")
                
                if not dry_run:
                    vasp_files = ['INCAR', 'INCAR-SC', 'INCAR-ELF', 'KPOINTS', 'POTCAR', 'POSCAR', 'CONTCAR',
                                  'OUTCAR', 'OSZICAR', 'OSZICAR-SC', 'vasprun.xml', 'vasprun.xml-SC',
                                  'EIGENVAL', 'DOSCAR', 'CHGCAR', 'CHGCAR-SC', 'CHG', 'WAVECAR', 'WAVECAR-SC',
                                  'WFULL', 'AECCAR0', 'AECCAR1', 'AECCAR2', 'PROCAR', 'LOCPOT', 'TMPCAR',
                                  'PCDAT', 'XDATCAR', 'REPORT', 'ELFCAR', 'PARCHG.tar.gz', 'IBZKPT',
                                  'SC_DONE', 'PARCHG_DONE', 'ELF_DONE', 'VASP_FAILED', 'VASP_DONE', 'job.sh',
                                  'generate_parchg_incars.py']
                    vasp_patterns = ['vasp_*.out', 'vasp_*.err', 'PARCHG-*', 'INCAR-PARCHG-*']
                    
                    removed_count = 0
                    for fname in vasp_files:
                        fpath = job_dir / fname
                        if fpath.exists():
                            fpath.unlink()
                            removed_count += 1
                    
                    for pattern in vasp_patterns:
                        for f in job_dir.glob(pattern):
                            f.unlink()
                            removed_count += 1
                    
                    print(f"      Removed {removed_count} files from SPE directory")
                    cleaned_dirs.append(str(job_dir))
            
            elif state in ['PARCHG_FAILED', 'ELF_FAILED']:
                # Since we reset to RELAX_DONE, clean entire SPE directory
                print(f"    - Removing all SPE directory contents (will redo entire SPE job)")
                relax_dir = Path(sdata['relax_dir'])
                print(f"    - Relax/CONTCAR will be preserved at: {relax_dir / 'CONTCAR'}")
                
                if not dry_run:
                    vasp_files = ['INCAR', 'INCAR-SC', 'INCAR-ELF', 'KPOINTS', 'POTCAR', 'POSCAR', 'CONTCAR',
                                  'OUTCAR', 'OSZICAR', 'OSZICAR-SC', 'vasprun.xml', 'vasprun.xml-SC',
                                  'EIGENVAL', 'DOSCAR', 'CHGCAR', 'CHGCAR-SC', 'CHG', 'WAVECAR', 'WAVECAR-SC',
                                  'WFULL', 'AECCAR0', 'AECCAR1', 'AECCAR2', 'PROCAR', 'LOCPOT', 'TMPCAR',
                                  'PCDAT', 'XDATCAR', 'REPORT', 'ELFCAR', 'PARCHG.tar.gz', 'IBZKPT',
                                  'SC_DONE', 'PARCHG_DONE', 'ELF_DONE', 'VASP_FAILED', 'VASP_DONE', 'job.sh',
                                  'generate_parchg_incars.py']
                    vasp_patterns = ['vasp_*.out', 'vasp_*.err', 'PARCHG-*', 'INCAR-PARCHG-*']
                    
                    removed_count = 0
                    for fname in vasp_files:
                        fpath = job_dir / fname
                        if fpath.exists():
                            fpath.unlink()
                            removed_count += 1
                    
                    for pattern in vasp_patterns:
                        for f in job_dir.glob(pattern):
                            f.unlink()
                            removed_count += 1
                    
                    print(f"      Removed {removed_count} files from SPE directory")
                    cleaned_dirs.append(str(job_dir))
    
    # Save updated database
    if not dry_run and reset_counts:
        backup_path = db_path.with_suffix('.json.bak')
        print(f"\nCreating backup: {backup_path}")
        shutil.copy2(db_path, backup_path)
        
        print(f"Saving updated database: {db_path}")
        with open(db_path, 'w') as f:
            json.dump(data, f, indent=2)
    
    # Print summary
    print("\n" + "="*80)
    print("Reset Summary")
    print("="*80)
    if reset_counts:
        for state, count in sorted(reset_counts.items()):
            new_state = reset_map[state][0]
            print(f"  {state} -> {new_state}: {count} structures")
        print(f"\nTotal reset: {sum(reset_counts.values())} structures")
    else:
        print("  No failed structures found matching criteria")
    
    if cleaned_dirs:
        print(f"\nCleaned directories: {len(cleaned_dirs)}")
    
    if dry_run:
        print("\nDRY RUN - No changes were made")
        print("Run without --dry-run to apply changes")
    else:
        print("\nChanges applied successfully!")
        print("Resume the workflow manager to retry failed jobs")
    
    print("="*80 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Reset failed VASP jobs for retry",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split('Usage:')[1]
    )
    parser.add_argument(
        '--db',
        type=str,
        default='./VASP_JOBS/workflow.json',
        help="Path to workflow.json database (default: ./VASP_JOBS/workflow.json)"
    )
    parser.add_argument(
        '--stage',
        type=str,
        choices=['RELAX', 'SC', 'PARCHG', 'ELF', 'relax', 'sc', 'parchg', 'elf'],
        help="Only reset specific stage (RELAX, SC, PARCHG, or ELF)"
    )
    parser.add_argument(
        '--clean',
        action='store_true',
        help="Remove VASP_FAILED markers from directories"
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help="Show what would be done without making changes"
    )
    parser.add_argument(
        '--list',
        action='store_true',
        help="List all failed structures and exit"
    )
    
    args = parser.parse_args()
    
    db_path = Path(args.db).expanduser()
    
    if not db_path.exists():
        print(f"Error: Database not found at {db_path}")
        return 1
    
    # List mode
    if args.list:
        with open(db_path, 'r') as f:
            data = json.load(f)
        list_failed_structures(data)
        return 0
    
    # Reset mode
    reset_failed_jobs(
        db_path=db_path,
        stage_filter=args.stage,
        clean=args.clean,
        dry_run=args.dry_run
    )
    
    return 0


if __name__ == '__main__':
    exit(main())

