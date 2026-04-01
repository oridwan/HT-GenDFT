#!/usr/bin/env python3
"""
Reset SC_FAILED, ELF_FAILED, PARCHG_FAILED jobs to RELAX_FAILED

This script resets all failed SPE calculation stages back to RELAX_FAILED,
allowing reset_failed_jobs.py to restart structures from scratch.

The SPE workflow is unified (SC → PARCHG → ELF in one job), so this script:
1. Creates timestamped backup of workflow.json
2. Changes SC_FAILED/ELF_FAILED/PARCHG_FAILED → RELAX_FAILED
3. Removes entire SPE directory (contains all stages: SC, PARCHG, ELF)
4. Creates/updates VASP_FAILED markers in RELAX directories
5. Clears spe_job_id and band gap fields (directory paths preserved for reuse)

After running this script, use reset_failed_jobs.py to retry from RELAX:
    python3 reset_failed_jobs.py --stage RELAX --clean > reset_log

Usage:
    python reset_failed_to_relax.py --dry-run    # Preview changes
    python reset_failed_to_relax.py              # Execute changes
"""

import json
import sys
import os
import shutil
import argparse
from datetime import datetime
from pathlib import Path


def load_workflow(workflow_path):
    """Load workflow.json database."""
    try:
        with open(workflow_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"ERROR: Could not load workflow database: {e}")
        sys.exit(1)


def save_workflow(workflow_path, data):
    """Save workflow.json database."""
    with open(workflow_path, 'w') as f:
        json.dump(data, f, indent=2)


def create_backup(workflow_path, dry_run=False):
    """Create timestamped backup of workflow.json."""
    workflow_dir = os.path.dirname(workflow_path)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = os.path.join(workflow_dir, f"workflow_{timestamp}.json.bak")
    
    if dry_run:
        print(f"  [DRY RUN] Would create: {backup_path}")
        return None
    
    shutil.copy2(workflow_path, backup_path)
    print(f"  Backup created: {backup_path}")
    
    # Verify backup
    try:
        with open(backup_path, 'r') as f:
            json.load(f)
        print(f"  Backup verified successfully")
        return backup_path
    except Exception as e:
        print(f"  ERROR: Backup verification failed: {e}")
        sys.exit(1)


def find_failed_structures(data):
    """Find all structures in failed states."""
    sc_failed = []
    elf_failed = []
    parchg_failed = []
    
    for struct_id, sdata in data['structures'].items():
        state = sdata.get('state')
        if state == 'SC_FAILED':
            sc_failed.append(struct_id)
        elif state == 'ELF_FAILED':
            elf_failed.append(struct_id)
        elif state == 'PARCHG_FAILED':
            parchg_failed.append(struct_id)
    
    return {
        'sc_failed': sorted(sc_failed),
        'elf_failed': sorted(elf_failed),
        'parchg_failed': sorted(parchg_failed),
        'all_failed': sorted(sc_failed + elf_failed + parchg_failed)
    }


def get_directory_size(path):
    """Get human-readable directory size."""
    try:
        total = 0
        for dirpath, dirnames, filenames in os.walk(path):
            for filename in filenames:
                filepath = os.path.join(dirpath, filename)
                if os.path.exists(filepath):
                    total += os.path.getsize(filepath)
        
        # Convert to human-readable
        for unit in ['B', 'KB', 'MB', 'GB']:
            if total < 1024.0:
                return f"{total:.1f}{unit}"
            total /= 1024.0
        return f"{total:.1f}TB"
    except:
        return "unknown"


def update_relax_marker(relax_dir, dry_run=False):
    """Update/create VASP_FAILED marker in RELAX directory."""
    if not relax_dir or not os.path.isdir(relax_dir):
        return False, f"RELAX directory not found: {relax_dir}"
    
    vasp_done = os.path.join(relax_dir, 'VASP_DONE')
    vasp_failed = os.path.join(relax_dir, 'VASP_FAILED')
    
    # If VASP_FAILED already exists, nothing to do
    if os.path.exists(vasp_failed):
        return False, "VASP_FAILED marker already exists"
    
    if os.path.exists(vasp_done):
        # Rename VASP_DONE → VASP_FAILED
        if dry_run:
            return True, "[DRY RUN] Would rename: VASP_DONE → VASP_FAILED"
        os.rename(vasp_done, vasp_failed)
        return True, "RELAX: VASP_DONE → VASP_FAILED"
    else:
        # No VASP_DONE found (likely timed out) - create VASP_FAILED marker
        if dry_run:
            return True, "[DRY RUN] Would create: VASP_FAILED marker (RELAX timed out)"
        with open(vasp_failed, 'w') as f:
            f.write("RELAX job failed or timed out - marked by reset_failed_to_relax.py\n")
        return True, "RELAX: Created VASP_FAILED marker (no VASP_DONE found)"


def remove_directory(dir_path, dir_type, dry_run=False):
    """Remove a directory (SC/ELF/PARCHG)."""
    if not dir_path or not os.path.isdir(dir_path):
        return False, None
    
    if dry_run:
        size = get_directory_size(dir_path)
        return True, f"[DRY RUN] Would remove {dir_type} directory ({size})"
    
    shutil.rmtree(dir_path)
    return True, f"Removed {dir_type} directory"


def process_structure(struct_id, sdata, dry_run=False):
    """Process a single failed structure."""
    results = {
        'relax_updated': False,
        'spe_removed': False,
        'messages': []
    }
    
    print(f"  Processing: {struct_id}")
    
    # Update RELAX marker
    relax_dir = sdata.get('relax_dir', '')
    updated, msg = update_relax_marker(relax_dir, dry_run)
    if updated:
        results['relax_updated'] = True
        results['messages'].append(f"    {msg}")
    else:
        results['messages'].append(f"    - RELAX: {msg}")
    
    # Remove unified SPE directory
    spe_dir = sdata.get('spe_dir', '')
    removed, msg = remove_directory(spe_dir, 'SPE', dry_run)
    if removed:
        results['spe_removed'] = True
        results['messages'].append(f"      {msg}")
    elif spe_dir:
        results['messages'].append(f"    - SPE directory not found: {spe_dir}")
    else:
        results['messages'].append(f"    - No SPE directory path in database")
    
    # Print messages
    for msg in results['messages']:
        print(msg)
    print()
    
    return results


def update_workflow_states(data, failed_structures, dry_run=False):
    """Update workflow.json to reset failed states to RELAX_FAILED."""
    if dry_run:
        return 0
    
    updated_count = 0
    failed_states = ['SC_FAILED', 'ELF_FAILED', 'PARCHG_FAILED']
    
    for struct_id in failed_structures:
        if struct_id in data['structures']:
            sdata = data['structures'][struct_id]
            
            if sdata.get('state') in failed_states:
                # Reset to RELAX_FAILED
                sdata['state'] = 'RELAX_FAILED'
                
                # Clear SPE job ID (unified job for SC-PARCHG-ELF)
                sdata['spe_job_id'] = None
                
                # Clear band gap info (will be recalculated in new SPE run)
                sdata['band_gap'] = None
                sdata['is_semiconductor'] = None
                
                updated_count += 1
    
    return updated_count


def verify_changes(workflow_path):
    """Verify the updated workflow.json."""
    data = load_workflow(workflow_path)
    
    state_counts = {}
    for sdata in data['structures'].values():
        state = sdata.get('state', 'UNKNOWN')
        state_counts[state] = state_counts.get(state, 0) + 1
    
    print("  Current state counts:")
    print(f"    ELF_FAILED: {state_counts.get('ELF_FAILED', 0)} (should be 0)")
    print(f"    SC_FAILED: {state_counts.get('SC_FAILED', 0)} (should be 0)")
    print(f"    PARCHG_FAILED: {state_counts.get('PARCHG_FAILED', 0)} (should be 0)")
    print(f"    RELAX_FAILED: {state_counts.get('RELAX_FAILED', 0)}")


def main():
    parser = argparse.ArgumentParser(
        description="Reset failed calculation stages to RELAX_FAILED",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        '--workflow',
        default='VASP_JOBS/workflow.json',
        help='Path to workflow.json (default: VASP_JOBS/workflow.json)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Preview changes without executing'
    )
    
    args = parser.parse_args()
    
    # Header
    if args.dry_run:
        print("=" * 67)
        print("DRY RUN MODE - No changes will be made")
        print("=" * 67)
        print()
    
    print("=" * 67)
    print("Reset Failed Jobs to RELAX_FAILED")
    print("=" * 67)
    print(f"Workflow database: {args.workflow}")
    print(f"Timestamp: {datetime.now().strftime('%Y%m%d_%H%M%S')}")
    print()
    print("Will reset: SC_FAILED, ELF_FAILED, PARCHG_FAILED → RELAX_FAILED")
    print()
    
    # Check workflow exists
    if not os.path.exists(args.workflow):
        print(f"ERROR: Workflow database not found: {args.workflow}")
        print(f"Current working directory: {os.getcwd()}")
        print(f"Make sure you're running from the correct directory")
        sys.exit(1)
    
    # Show working directory for troubleshooting
    print(f"Working directory: {os.getcwd()}")
    print()
    
    # Step 1: Create backup
    print("Step 1: Creating backup...")
    backup_path = create_backup(args.workflow, args.dry_run)
    print()
    
    # Load workflow data
    data = load_workflow(args.workflow)
    
    # Step 2: Find failed structures
    print("Step 2: Finding failed structures...")
    failed = find_failed_structures(data)
    total_count = len(failed['all_failed'])
    
    if total_count == 0:
        print("  No failed structures found.")
        sys.exit(0)
    
    print(f"  Found {total_count} failed structures:")
    print(f"    SC_FAILED: {len(failed['sc_failed'])}")
    print(f"    ELF_FAILED: {len(failed['elf_failed'])}")
    print(f"    PARCHG_FAILED: {len(failed['parchg_failed'])}")
    print()
    
    # Step 3: Process each structure
    print("Step 3: Processing structures...")
    total_relax_updated = 0
    total_spe_removed = 0
    
    for struct_id in failed['all_failed']:
        sdata = data['structures'].get(struct_id, {})
        results = process_structure(struct_id, sdata, args.dry_run)
        
        if results['relax_updated']:
            total_relax_updated += 1
        if results['spe_removed']:
            total_spe_removed += 1
    
    # Step 4: Update workflow.json
    print("Step 4: Updating workflow.json...")
    if args.dry_run:
        print(f"  [DRY RUN] Would update {total_count} structures → RELAX_FAILED")
        print(f"  [DRY RUN] Would clear job_id fields (keeping directory paths)")
    else:
        updated_count = update_workflow_states(data, failed['all_failed'], args.dry_run)
        save_workflow(args.workflow, data)
        
        # Verify saved JSON is valid
        try:
            verify_data = load_workflow(args.workflow)
            print(f"    Updated {updated_count} structures → RELAX_FAILED")
            print(f"    workflow.json updated successfully")
        except Exception as e:
            print(f"  ERROR: Updated workflow.json is invalid: {e}")
            if backup_path:
                print(f"  Restoring from backup...")
                shutil.copy2(backup_path, args.workflow)
            sys.exit(1)
    print()
    
    # Step 5: Verify changes
    if not args.dry_run:
        print("Step 5: Verifying changes...")
        verify_changes(args.workflow)
        print()
    
    # Summary
    print("=" * 67)
    print("Summary")
    print("=" * 67)
    
    if args.dry_run:
        print("[DRY RUN] No changes were made.")
        print()
        print(f"Would process: {total_count} structures")
        print(f"Would update: RELAX directories (VASP_DONE → VASP_FAILED)")
        print(f"Would remove: SPE directories (unified SC-PARCHG-ELF)")
        print(f"Would update workflow.json: {total_count} structures → RELAX_FAILED")
        print()
        print("To execute these changes, run without --dry-run:")
        print("  python reset_failed_to_relax.py")
    else:
        print(f"  Backup created: {backup_path}")
        print(f"  Processed structures: {total_count}")
        print(f"  RELAX directories updated: {total_relax_updated}")
        print(f"  SPE directories removed: {total_spe_removed}")
        print(f"  workflow.json updated: SC_FAILED/ELF_FAILED/PARCHG_FAILED → RELAX_FAILED")
        print()
        print("All operations completed successfully!")
        print()
        print("Next step: Run reset_failed_jobs.py to retry from RELAX")
        print("  python3 reset_failed_jobs.py --stage RELAX --clean > reset_log")
        print()
        print("To restore from backup if needed:")
        print(f"  cp {backup_path} {args.workflow}")
    
    print("=" * 67)


if __name__ == '__main__':
    main()

