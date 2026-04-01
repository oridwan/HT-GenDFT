#!/usr/bin/env bash
set -euo pipefail

TARGET_USER="${1:-$USER}"

if ! states="$(squeue -h -u "$TARGET_USER" -o "%t" 2>/dev/null)"; then
    echo "Error: failed to query SLURM with squeue for user '$TARGET_USER'" >&2
    exit 1
fi

running="$(awk '$1=="R"{c++} END{print c+0}' <<< "$states")"
pending="$(awk '$1=="PD"{c++} END{print c+0}' <<< "$states")"
total="$(awk 'NF{c++} END{print c+0}' <<< "$states")"

echo "User: $TARGET_USER"
echo "Total jobs: $total"
echo "Running: $running"
echo "Pending: $pending"

if (( pending > 0 )); then
    echo ""
    echo "Pending reasons:"
    if ! reasons="$(squeue -h -u "$TARGET_USER" -t PD -o "%r" 2>/dev/null)"; then
        echo "  (Could not read pending reasons)"
    else
        printf '%s\n' "$reasons" \
            | sort \
            | uniq -c \
            | sort -nr \
            | awk '{printf "  %4d  %s\n", $1, substr($0, index($0, $2))}'
    fi
fi
