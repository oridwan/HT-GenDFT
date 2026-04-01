#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash git_push_scripts_only.sh [options] "Commit message"

Options:
  --no-push         Create the commit but do not push it.
  --dry-run         Show which script files would be committed.
  --remote NAME     Remote to push to. Default: origin
  --branch NAME     Branch to push to. Default: current branch
  -h, --help        Show this help.

Behavior:
  - Includes only changed .py and .sh files.
  - Includes tracked deletions of .py and .sh files.
  - Leaves unrelated staged changes untouched.
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

remote="origin"
branch=""
push_after_commit=1
dry_run=0

while (($#)); do
  case "$1" in
    --no-push)
      push_after_commit=0
      shift
      ;;
    --dry-run)
      dry_run=1
      shift
      ;;
    --remote)
      [[ $# -ge 2 ]] || die "--remote requires a value"
      remote="$2"
      shift 2
      ;;
    --branch)
      [[ $# -ge 2 ]] || die "--branch requires a value"
      branch="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      break
      ;;
    -*)
      die "unknown option: $1"
      ;;
    *)
      break
      ;;
  esac
done

[[ $# -ge 1 ]] || die "commit message is required"
commit_message="$1"

git rev-parse --is-inside-work-tree >/dev/null 2>&1 || die "not inside a git repository"

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

index_lock="$(git rev-parse --git-path index.lock)"
if [[ -e "$index_lock" ]]; then
  die "git index is locked at $index_lock; remove it after confirming no other git process is using this repository"
fi

if [[ -z "$branch" ]]; then
  branch="$(git branch --show-current)"
fi
[[ -n "$branch" ]] || die "detached HEAD; pass --branch explicitly"

declare -a raw_paths=()
declare -A seen=()
declare -a script_paths=()

while IFS= read -r -d '' entry; do
  status="${entry:0:2}"
  path="${entry:3}"
  raw_paths+=("$path")

  if [[ "$status" == *R* || "$status" == *C* ]]; then
    IFS= read -r -d '' path || true
    [[ -n "$path" ]] && raw_paths+=("$path")
  fi
done < <(git status --porcelain=v1 -z)

for path in "${raw_paths[@]}"; do
  [[ -n "$path" ]] || continue
  [[ "$path" == *.py || "$path" == *.sh ]] || continue
  if [[ -z "${seen[$path]+x}" ]]; then
    seen["$path"]=1
    script_paths+=("$path")
  fi
done

if [[ ${#script_paths[@]} -eq 0 ]]; then
  echo "No changed .py or .sh files found."
  exit 0
fi

tmp_index="$(mktemp "${TMPDIR:-/tmp}/git-script-only-index.XXXXXX")"
cleanup() {
  rm -f "$tmp_index"
}
trap cleanup EXIT

GIT_INDEX_FILE="$tmp_index" git read-tree HEAD
GIT_INDEX_FILE="$tmp_index" git add -A -- "${script_paths[@]}"

if GIT_INDEX_FILE="$tmp_index" git diff --cached --quiet; then
  echo "No script-only diff relative to HEAD."
  exit 0
fi

echo "Script paths selected:"
printf '  %s\n' "${script_paths[@]}"

if ((dry_run)); then
  exit 0
fi

parent_commit="$(git rev-parse HEAD)"
tree_id="$(GIT_INDEX_FILE="$tmp_index" git write-tree)"
commit_id="$(printf '%s\n' "$commit_message" | git commit-tree "$tree_id" -p "$parent_commit")"

git update-ref "refs/heads/$branch" "$commit_id" "$parent_commit"
git reset --quiet HEAD -- "${script_paths[@]}"

echo "Created commit $commit_id on branch $branch"

if ((push_after_commit)); then
  git push "$remote" "$branch"
fi
