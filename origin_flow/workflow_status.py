#!/usr/bin/env python3
"""
Workflow Status Checker - Query the workflow database without starting jobs
"""

import sys
import json
import argparse
from pathlib import Path
from datetime import datetime
from collections import defaultdict


def load_db(db_path):
    """Load workflow database."""
    db_path = Path(db_path)
    if not db_path.exists():
        print(f"Error: Database not found: {db_path}")
        sys.exit(1)
    
    with open(db_path, 'r') as f:
        return json.load(f)


def print_summary(data):
    """Print overall workflow summary."""
    structures = data['structures']
    config = data.get('config', {})
    
    print("\n" + "="*70)
    print("WORKFLOW SUMMARY")
    print("="*70)
    
    if config:
        print(f"Results dir: {config.get('results_dir', 'N/A')}")
        print(f"Output dir: {config.get('output_dir', 'N/A')}")
        print(f"Max concurrent: {config.get('max_concurrent', 'N/A')}")
    
    print(f"\nTotal structures: {len(structures)}")
    
    # Count by state
    state_counts = defaultdict(int)
    for s in structures.values():
        state_counts[s['state']] += 1
    
    print("\nStatus breakdown:")
    for state in sorted(state_counts.keys()):
        count = state_counts[state]
        pct = 100 * count / len(structures) if structures else 0
        print(f"  {state:20s}: {count:4d} ({pct:5.1f}%)")
    
    # Running count
    running_states = ['RELAX_RUNNING', 'SC_RUNNING', 'PARCHG_RUNNING', 'ELF_RUNNING']
    running = sum(state_counts[s] for s in running_states)
    print(f"\n  Currently running: {running}")
    
    # Completion stats
    completed = state_counts['ELF_DONE']
    failed = sum(state_counts[s] for s in state_counts if 'FAILED' in s)
    pending = state_counts['PENDING'] + sum(state_counts[s] for s in state_counts if 'DONE' in s and s != 'ELF_DONE')
    
    print(f"\nProgress:")
    print(f"  Completed: {completed}/{len(structures)} ({100*completed/len(structures) if structures else 0:.1f}%)")
    print(f"  Failed: {failed}")
    print(f"  Pending/In-progress: {pending + running}")
    
    print("="*70 + "\n")


def print_details(data, state_filter=None, composition_filter=None, limit=None):
    """Print detailed structure information."""
    structures = data['structures']
    
    # Filter structures
    filtered = []
    for sid, s in structures.items():
        if state_filter and s['state'] != state_filter:
            continue
        if composition_filter and s['composition'] != composition_filter:
            continue
        filtered.append((sid, s))
    
    if not filtered:
        print("No structures match the filters.")
        return
    
    # Sort by last updated
    filtered.sort(key=lambda x: x[1]['last_updated'], reverse=True)
    
    if limit:
        filtered = filtered[:limit]
    
    print(f"\n{'ID':<20s} {'State':<20s} {'Composition':<15s} {'Last Updated':<20s}")
    print("-"*75)
    
    for sid, s in filtered:
        last_update = s['last_updated'][:19] if s['last_updated'] else 'N/A'
        print(f"{sid:<20s} {s['state']:<20s} {s['composition']:<15s} {last_update:<20s}")
        
        if s.get('error'):
            print(f"  Error: {s['error']}")
    
    print()


def print_compositions(data):
    """Print composition-level summary."""
    structures = data['structures']
    
    comp_stats = defaultdict(lambda: defaultdict(int))
    for s in structures.values():
        comp = s['composition']
        state = s['state']
        comp_stats[comp][state] += 1
    
    print("\n" + "="*70)
    print("COMPOSITION SUMMARY")
    print("="*70)
    print(f"\n{'Composition':<15s} {'Total':>6s} {'Done':>6s} {'Failed':>6s} {'Running':>6s} {'Pending':>6s}")
    print("-"*70)
    
    for comp in sorted(comp_stats.keys()):
        stats = comp_stats[comp]
        total = sum(stats.values())
        done = stats.get('ELF_DONE', 0)
        failed = sum(stats[s] for s in stats if 'FAILED' in s)
        running = sum(stats[s] for s in ['RELAX_RUNNING', 'SC_RUNNING', 'PARCHG_RUNNING', 'ELF_RUNNING'])
        pending = total - done - failed - running
        
        print(f"{comp:<15s} {total:6d} {done:6d} {failed:6d} {running:6d} {pending:6d}")
    
    print("="*70 + "\n")


def print_failed(data):
    """Print all failed structures with details."""
    structures = data['structures']
    
    failed = [(sid, s) for sid, s in structures.items() 
              if 'FAILED' in s['state']]
    
    if not failed:
        print("\nNo failed structures.")
        return
    
    print("\n" + "="*70)
    print("FAILED STRUCTURES")
    print("="*70)
    
    for sid, s in sorted(failed):
        print(f"\n{sid}")
        print(f"  State: {s['state']}")
        print(f"  Composition: {s['composition']}")
        print(f"  Last updated: {s['last_updated'][:19]}")
        
        if s.get('error'):
            print(f"  Error: {s['error']}")
        
        # Print job directories
        if s['state'] == 'RELAX_FAILED':
            print(f"  Check: {s['relax_dir']}")
        elif s['state'] in ['SC_FAILED', 'PARCHG_FAILED', 'ELF_FAILED']:
            print(f"  Check: {s['spe_dir']}")
    
    print("="*70 + "\n")


def print_running(data):
    """Print currently running jobs."""
    structures = data['structures']
    
    running = [(sid, s) for sid, s in structures.items() 
               if 'RUNNING' in s['state']]
    
    if not running:
        print("\nNo running jobs.")
        return
    
    print("\n" + "="*70)
    print("RUNNING JOBS")
    print("="*70)
    print(f"\n{'ID':<20s} {'Stage':<15s} {'Job ID':<15s} {'Composition':<15s}")
    print("-"*70)
    
    for sid, s in sorted(running):
        state = s['state']
        
        if state == 'RELAX_RUNNING':
            stage = 'Relax'
            job_id = s.get('relax_job_id', 'N/A')
        elif state in ['SC_RUNNING', 'PARCHG_RUNNING', 'ELF_RUNNING']:
            # All SPE stages use the same job ID
            if state == 'SC_RUNNING':
                stage = 'SPE (SC)'
            elif state == 'PARCHG_RUNNING':
                stage = 'SPE (PARCHG)'
            elif state == 'ELF_RUNNING':
                stage = 'SPE (ELF)'
            job_id = s.get('spe_job_id', 'N/A')
        else:
            stage = 'Unknown'
            job_id = 'N/A'
        
        print(f"{sid:<20s} {stage:<15s} {job_id:<15s} {s['composition']:<15s}")
    
    print("="*70 + "\n")


def print_completed(data, limit=10):
    """Print recently completed structures."""
    structures = data['structures']
    
    completed = [(sid, s) for sid, s in structures.items() 
                 if s['state'] == 'ELF_DONE']
    
    if not completed:
        print("\nNo completed structures.")
        return
    
    # Sort by last updated
    completed.sort(key=lambda x: x[1]['last_updated'], reverse=True)
    
    print("\n" + "="*70)
    print(f"COMPLETED STRUCTURES (showing {min(limit, len(completed))} most recent)")
    print("="*70)
    print(f"\n{'ID':<20s} {'Composition':<15s} {'Completed':<20s}")
    print("-"*70)
    
    for sid, s in completed[:limit]:
        last_update = s['last_updated'][:19]
        print(f"{sid:<20s} {s['composition']:<15s} {last_update:<20s}")
    
    print("="*70 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Check VASP workflow status"
    )
    parser.add_argument(
        'db_path',
        nargs='?',
        default='workflow.json',
        help="Path to workflow database (default: workflow.json)"
    )
    parser.add_argument(
        '--details',
        action='store_true',
        help="Show detailed structure list"
    )
    parser.add_argument(
        '--state',
        type=str,
        help="Filter by state (e.g., PENDING, RELAX_RUNNING, ELF_DONE)"
    )
    parser.add_argument(
        '--composition',
        type=str,
        help="Filter by composition"
    )
    parser.add_argument(
        '--failed',
        action='store_true',
        help="Show failed structures"
    )
    parser.add_argument(
        '--running',
        action='store_true',
        help="Show running jobs"
    )
    parser.add_argument(
        '--completed',
        action='store_true',
        help="Show completed structures"
    )
    parser.add_argument(
        '--compositions',
        action='store_true',
        help="Show composition-level summary"
    )
    parser.add_argument(
        '--limit',
        type=int,
        default=20,
        help="Limit number of results shown (default: 20)"
    )
    
    args = parser.parse_args()
    
    # Load database
    data = load_db(args.db_path)
    
    # Print summary
    print_summary(data)
    
    # Additional views
    if args.compositions:
        print_compositions(data)
    
    if args.failed:
        print_failed(data)
    
    if args.running:
        print_running(data)
    
    if args.completed:
        print_completed(data, args.limit)
    
    if args.details:
        print_details(
            data, 
            state_filter=args.state,
            composition_filter=args.composition,
            limit=args.limit
        )


if __name__ == '__main__':
    main()

