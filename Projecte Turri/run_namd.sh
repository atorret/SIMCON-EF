#!/bin/bash
# Runs the CHARMM-GUI NAMD steps of ONE namd folder, one after another:
#   6.1 6.2 6.3 6.4 6.5 6.6 (equilibration)  ->  7 (production, log in output.log)
#
# Usage (in WSL):
#   bash run_namd.sh <namd-folder>              # runs every step that is not finished yet
#   bash run_namd.sh <namd-folder> 6.6          # (re)runs from 6.6 to 7, even if finished
#   bash run_namd.sh <namd-folder> 6.6 6.6      # (re)runs only 6.6
#
# Without a first step it RESUMES: steps whose log ends correctly are skipped,
# and everything from the first unfinished step onwards is run.
#
# The run happens in a copy on the WSL Linux disk (~/namd-runs/...), because writing
# to the Windows disk (/mnt/c) can fail with "Permission denied" while Windows is
# reading the files. Results are copied back to <namd-folder> after every step.

NAMD="${NAMD:-$HOME/namd/namd3}"            # Linux NAMD 3 (CUDA) installed in WSL
NAMD_OPTS="+p10 +devices 0"        # 10 CPU threads + GPU (leaves cores free for Windows)

STEPS=(6.1 6.2 6.3 6.4 6.5 6.6 7)

log() { echo "[$(date '+%F %T')] $*"; }

[ -n "$1" ] || { echo "Usage: bash run_namd.sh <namd-folder> [first-step] [last-step]"; exit 1; }
SRC="$(cd "$1" 2>/dev/null && pwd)" || { log "ERROR: folder not found: $1"; exit 1; }
FIRST="$2"
LAST="${3:-7}"
[ -x "$NAMD" ] || { log "ERROR: NAMD not found at $NAMD"; exit 1; }

# Unique work folder name, e.g. charmm-ORG_charmm-gui-9562151784_namd-TEMP-283
PARENT="$(dirname "$SRC")"
TAG="$(basename "$(dirname "$PARENT")")_$(basename "$PARENT")_$(basename "$SRC")"
WORK_PARENT="$HOME/namd-runs/$TAG"
WORK="$WORK_PARENT/$(basename "$SRC")"

input_for() { [ "$1" = "7" ] && echo "step7_production" || echo "step$1_equilibration"; }
log_for()   { [ "$1" = "7" ] && echo "output.log"       || echo "step$1_equilibration.out"; }

# NAMD prints "End of program" even after a crash; "WallClock:" only on a normal finish.
finished() { [ -f "$1" ] && grep -q "^WallClock:" "$1" && ! grep -q "FATAL ERROR" "$1"; }

copy_back() {
    rsync -a "$WORK"/ "$SRC"/ \
        && log "Results copied back to $SRC" \
        || log "WARNING: could not copy results back. They are in $WORK"
}

log "===== $SRC"
log "Fresh copy -> $WORK"
rm -rf "$WORK_PARENT"
mkdir -p "$WORK_PARENT"
rsync -a --exclude 'run_all.log' "$SRC"/ "$WORK"/ || { log "ERROR: copy failed"; exit 1; }
# step6.1 reads ../step5_assembly.str from the CHARMM-GUI folder
[ -f "$PARENT/step5_assembly.str" ] && cp "$PARENT/step5_assembly.str" "$WORK_PARENT"/
cd "$WORK" || exit 1

running=0
[ -z "$FIRST" ] && resume=1 || resume=0
for s in "${STEPS[@]}"; do
    name=$(input_for "$s")
    out=$(log_for "$s")

    if [ $running -eq 0 ]; then
        if [ $resume -eq 1 ]; then
            if finished "$out"; then log "Step $s already finished, skipping"; continue; fi
            running=1
        else
            [ "$s" = "$FIRST" ] && running=1 || continue
        fi
    fi

    log "Starting step $s  ($name.inp -> $out)"
    "$NAMD" $NAMD_OPTS "$name.inp" > "$out" 2>&1

    if ! finished "$out"; then
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

log "Done: $SRC"
