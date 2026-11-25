#!/bin/bash
# Runs the NAMD steps one after another (equilibration 6.1-6.6, then production 7).
# Each step starts only if the previous one finished correctly.
#
# The simulation runs in a copy of this folder on the WSL Linux disk
# (~/namd-runs/<folder name>), because writing to the Windows disk (/mnt/c)
# can fail with "Permission denied" when Windows (VS Code, git, antivirus)
# is reading the same files. Results are copied back here after every step.
#
# Usage (from WSL, inside this folder):
#   bash run_all.sh            # runs 6.1 6.2 6.3 6.4 6.5 6.6 7
#   bash run_all.sh 6.4        # starts at 6.4 (6.1-6.3 already done)
#   bash run_all.sh 6.4 6.6    # runs only 6.4 to 6.6
#
# To keep it running after closing the terminal:
#   nohup bash run_all.sh > run_all.log 2>&1 &

SRC="$(cd "$(dirname "$0")" && pwd)"
WORK="$HOME/namd-runs/$(basename "$SRC")"

# Linux NAMD 3 installed in WSL (~/namd). The ./namd2 in this folder is a macOS binary.
NAMD="$HOME/namd/namd3"
NAMD_OPTS="+p10 +devices 0"     # 10 CPU threads + GPU (leaves cores free for Windows)

STEPS=(6.1 6.2 6.3 6.4 6.5 6.6 7)
FIRST="${1:-6.1}"
LAST="${2:-7}"

log() { echo "[$(date '+%F %T')] $*"; }

if [ ! -x "$NAMD" ]; then
    log "ERROR: NAMD not found at $NAMD"; exit 1
fi

copy_back() {
    rsync -a --exclude 'run_all.log' "$WORK"/ "$SRC"/ \
        && log "Results copied back to $SRC" \
        || log "WARNING: could not copy results back. They are in $WORK"
}

log "Copying $SRC -> $WORK"
mkdir -p "$WORK"
rsync -a --delete --exclude 'run_all.log' "$SRC"/ "$WORK"/ || { log "ERROR: copy failed"; exit 1; }
cd "$WORK" || exit 1

input_for() {   # step number -> input file name
    if [ "$1" = "7" ]; then echo "step7_production"; else echo "step$1_equilibration"; fi
}

running=0
for s in "${STEPS[@]}"; do
    [ "$s" = "$FIRST" ] && running=1
    [ $running -eq 0 ] && continue

    name=$(input_for "$s")
    out="$name.out"
    [ "$s" = "7" ] && out="output.log"   # production log name used in the project

    log "Starting step $s  ($name.inp -> $out)"
    $NAMD $NAMD_OPTS "$name.inp" > "$out" 2>&1

    # NAMD prints "End of program" even after a crash; "WallClock:" only on a normal finish.
    if grep -q "FATAL ERROR" "$out" || ! grep -q "^WallClock:" "$out"; then
        log "ERROR: step $s did not finish correctly. Check $out"
        grep -E "ERROR" "$out" | head -5
        tail -n 5 "$out"
        copy_back
        exit 1
    fi
    log "Step $s finished. $(grep '^WallClock:' "$out")"
    copy_back

    [ "$s" = "$LAST" ] && break
done

log "All requested steps done."
