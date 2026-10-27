#!/usr/bin/env bash
# Three-way merge dry run of the qe-pert (QE 7.3) DFPT+SOC+U changes onto
# this QE 7.6 tree. Nothing in the tree is modified; results go to OUT_DIR.
#   base   = QE 7.3 file at the pristine qe-pert commit
#   ours   = qe-pert file at SRC_REF
#   theirs = this tree's file at TARGET_REF (default: pristine 7.6 baseline),
#            following the known 7.3 -> 7.6 relocations below
set -euo pipefail

usage() {
  cat >&2 <<EOF
Usage: bash $0 OUT_DIR [--pert DIR] [--base REF] [--src REF] [--target REF|WORKTREE]
Defaults: --pert ~/bin/qe-pert --base 347e370 --src ph-acc --target 5e4b172
EOF
  exit 2
}

[[ $# -ge 1 ]] || usage
out=$1; shift
pert=$HOME/bin/qe-pert; base_ref=347e370; src_ref=ph-acc; target_ref=5e4b172
while [[ $# -gt 0 ]]; do
  case $1 in
    --pert) pert=$2; shift 2 ;;
    --base) base_ref=$2; shift 2 ;;
    --src) src_ref=$2; shift 2 ;;
    --target) target_ref=$2; shift 2 ;;
    *) usage ;;
  esac
done
tree=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)

# 7.3 path -> 7.6 path. An empty target means the 7.3 file has no 7.6 file.
declare -A moved=(
  [PHonon/PH/dnsq_scf.f90]=LR_Modules/dnsq_scf.f90
  [PHonon/PH/sym_dns.f90]=LR_Modules/sym_dns.f90
  [PHonon/PH/syme_dns.f90]=
)

target_file() {  # $1 = 7.6 path; prints content or nothing
  if [[ $target_ref == WORKTREE ]]; then
    [[ -f $tree/$1 ]] && cat "$tree/$1" || true
  else
    git -C "$tree" show "$target_ref:$1" 2>/dev/null || true
  fi
}

mkdir -p "$out/files"
{
  echo "# pert=$pert base=$base_ref src=$src_ref($(git -C "$pert" rev-parse --short "$src_ref")) target=$target_ref"
  printf '%-36s %-30s %6s %6s %5s\n' FILE_73 FILE_76 dPERT d73_76 CONFL
} > "$out/summary.txt"

git -C "$pert" diff --name-only "$base_ref" "$src_ref" -- LR_Modules PHonon/PH PW/src Modules upflib |
  grep -E '\.(f90|F90|c|h)$|/Makefile$' |
while read -r f73; do
  if [[ -v moved[$f73] ]]; then f76=${moved[$f73]}; else f76=$f73; fi
  stem=$out/files/${f73//\//_}
  git -C "$pert" show "$base_ref:$f73" > "$stem.base" 2>/dev/null || : > "$stem.base"
  git -C "$pert" show "$src_ref:$f73" > "$stem.ours"
  if [[ -n $f76 ]]; then target_file "$f76" > "$stem.theirs"; else : > "$stem.theirs"; fi
  dp=$( (diff "$stem.base" "$stem.ours" || true) | grep -c '^[<>]' || true)
  du=$( (diff "$stem.base" "$stem.theirs" || true) | grep -c '^[<>]' || true)
  cp "$stem.ours" "$stem.merged"
  set +e
  git merge-file -q -L qe-pert -L qe-7.3 -L qe-7.6 "$stem.merged" "$stem.base" "$stem.theirs"
  nconf=$?
  set -e
  printf '%-36s %-30s %6s %6s %5s\n' "$f73" "${f76:-(none)}" "$dp" "$du" "$nconf" >> "$out/summary.txt"
done
cat "$out/summary.txt"
echo "Per-file base/ours/theirs/merged copies: $out/files" >&2
