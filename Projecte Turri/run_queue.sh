#!/bin/bash
# Runs every simulation of the project one after another:
#   1. N-64   (303.15 K, steps 6.1 -> 7)
#   2. ORG    (128 lipids, 303.15 K, steps 6.1 -> 7)
#   3. N-192  (303.15 K, steps 6.1 -> 7)
#   4. ORG at 283 K  (copy of ORG/namd -> namd-TEMP-283, steps 6.6 and 7 only)
#   5. ORG at 323 K  (copy of ORG/namd -> namd-TEMP-323, steps 6.6 and 7 only)
#
# Usage (in WSL, from the "Projecte Turri" folder):
#   nohup bash run_queue.sh > run_queue.log 2>&1 &
#
# If it is interrupted, run the same command again: finished steps are skipped
# and it continues from the first unfinished step.
# If one simulation fails, the queue stops there (the next ones need a correct result).

cd "$(dirname "$0")" || exit 1
ROOT="$(pwd)"

N64="charmm-N-64/charmm-gui-9562164596"
ORG="charmm-ORG/charmm-gui-9562151784"
N192="charmm-N-192/charmm-gui-9562196952"
TEMPS=(283 323)

log() { echo "[$(date '+%F %T')] $*"; }
run() { bash "$ROOT/run_namd.sh" "$@" || { log "QUEUE STOPPED: $1 failed"; exit 1; }; }

# Copies ORG/namd (with 6.1-6.5 finished at 303.15 K) to namd-TEMP-X and sets
# temperature X in steps 6.6 and 7. Done only once: an existing copy is kept,
# so an interrupted temperature run resumes instead of starting over.
make_temp() {
    local T=$1 dst="$ORG/namd-TEMP-$1"
    if [ -d "$dst" ]; then log "$dst already exists, reusing it"; return; fi
    log "Creating $dst"
    rsync -a \
        --exclude 'step6.6_equilibration.*' --exclude 'step7_production.*' \
        --exclude 'output.log' --exclude 'run_all.log' \
        "$ORG/namd"/ "$dst"/ || { log "ERROR: copy failed"; exit 1; }
    cp "$ORG/namd/step6.6_equilibration.inp" "$ORG/namd/step7_production.inp" "$dst"/
    sed -i -E "s/^set temp +[0-9.]+;/set temp                $T;/" \
        "$dst/step6.6_equilibration.inp" "$dst/step7_production.inp"
    grep -H "^set temp" "$dst/step6.6_equilibration.inp" "$dst/step7_production.inp"
}

log "QUEUE START"
run "$N64/namd"
run "$ORG/namd"
run "$N192/namd"
for T in "${TEMPS[@]}"; do
    make_temp "$T"
    run "$ORG/namd-TEMP-$T"      # resumes: 6.1-6.5 are already finished, so it runs 6.6 and 7
done
log "QUEUE FINISHED: all simulations done"
