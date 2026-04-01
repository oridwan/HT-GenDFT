#!/usr/bin/env bash
set -euo pipefail

# cleanup_structures_by_date.sh
# List and optionally remove directories matching a pattern (default '*_structures')
# within a specified date range. Creates a Markdown summary file when listing.
#
# Usage examples:
#   # Dry-run: list dirs modified Feb 4–5, 2026
#   ./cleanup_structures_by_date.sh --start-date 2026-02-04 --end-date 2026-02-05
#
#   # Delete those directories after confirmation
#   ./cleanup_structures_by_date.sh --start-date 2026-02-04 --end-date 2026-02-05 --delete
#
#   # Non-interactive delete (use carefully)
#   ./cleanup_structures_by_date.sh --start-date 2026-02-04 --end-date 2026-02-05 --delete --yes
#

DIR="$(pwd)/results/quaternary_hydride"
PATTERN='*_structures'
START_DATE=""
END_DATE=""
DRY_RUN=1
DO_DELETE=0
ASSUME_YES=0
OUTPUT_FILE=""

print_usage() {
    cat <<EOF
Usage: $0 --start-date YYYY-MM-DD --end-date YYYY-MM-DD [options]

Options:
  --dir PATH          Directory to search (default: ./results/quaternary_hydride)
  --pattern GLOB      Filename pattern (default: '*_structures')
  --start-date DATE   Start date (inclusive) in YYYY-MM-DD
  --end-date DATE     End date (inclusive) in YYYY-MM-DD
  --list-only         (default) only list and write summary
  --delete            Delete matched directories (requires confirmation)
  --yes               Skip confirmation when deleting (dangerous)
  --output FILE       Save listing to FILE (default: modified_dirs_<start>_to_<end>.md in target dir)
  -h, --help          Show this help and exit

Examples:
  $0 --start-date 2026-02-04 --end-date 2026-02-05
  $0 --start-date 2026-02-04 --end-date 2026-02-05 --delete
  $0 --start-date 2026-02-04 --end-date 2026-02-05 --delete --yes
EOF
}

if [ $# -eq 0 ]; then
    print_usage
    exit 1
fi

while [ $# -gt 0 ]; do
    case "$1" in
        --dir) DIR="$2"; shift 2;;
        --pattern) PATTERN="$2"; shift 2;;
        --start-date) START_DATE="$2"; shift 2;;
        --end-date) END_DATE="$2"; shift 2;;
        --list-only) DRY_RUN=1; shift 1;;
        --delete) DO_DELETE=1; DRY_RUN=0; shift 1;;
        --yes) ASSUME_YES=1; shift 1;;
        --output) OUTPUT_FILE="$2"; shift 2;;
        -h|--help) print_usage; exit 0;;
        *) echo "Unknown option: $1"; print_usage; exit 2;;
    esac
done

if [ -z "$START_DATE" ] || [ -z "$END_DATE" ]; then
    echo "ERROR: --start-date and --end-date are required"
    print_usage
    exit 2
fi

# Ensure target directory exists
if [ ! -d "$DIR" ]; then
    echo "ERROR: Directory not found: $DIR"
    exit 2
fi

# Default output file if not supplied
if [ -z "$OUTPUT_FILE" ]; then
    OUTPUT_FILE="$DIR/modified_dirs_${START_DATE}_to_${END_DATE}.md"
fi

echo "Searching in: $DIR"
echo "Pattern: $PATTERN"
echo "Date range: $START_DATE (inclusive) -> $END_DATE (inclusive)"
echo "Output: $OUTPUT_FILE"

# Build and run the find command. Use end-of-day on END_DATE to include the whole day.
FIND_CMD=(find "$DIR" -maxdepth 1 -type d -name "$PATTERN" \
    -newermt "$START_DATE 00:00:00" ! -newermt "$END_DATE 23:59:59")

echo "\nMatched directories:"
"${FIND_CMD[@]}" -printf '%TY-%Tm-%Td %TT %p\n'

# Save listing to markdown file
{
    echo "# Modified directories ($START_DATE to $END_DATE)"
    echo "Generated on: $(date '+%F %T %z')"
    echo
    "${FIND_CMD[@]}" -printf ' - %TY-%Tm-%Td %TT %p\n'
} > "$OUTPUT_FILE"

echo "\nSaved listing to: $OUTPUT_FILE"

if [ "$DO_DELETE" -eq 1 ]; then
    # Count matches
    MATCH_COUNT=$("${FIND_CMD[@]}" | wc -l)
    echo "\nDelete requested. $MATCH_COUNT directories matched."

    if [ $MATCH_COUNT -eq 0 ]; then
        echo "Nothing to delete. Exiting."
        exit 0
    fi

    if [ $ASSUME_YES -ne 1 ]; then
        read -r -p "Proceed to delete these $MATCH_COUNT directories? (y/N) " CONFIRM
        case "$CONFIRM" in
            [yY]|[yY][eE][sS]) ;;
            *) echo "Aborted by user."; exit 0;;
        esac
    else
        echo "--yes given: proceeding without confirmation"
    fi

    # Perform deletion
    echo "Deleting..."
    # Use -exec rm -rf {} + to batch deletes safely
    "${FIND_CMD[@]}" -exec rm -rf {} +
    echo "Deletion completed."
fi

echo "Done."
