# How to run the simulations

All simulations run in **WSL (Ubuntu)** with NAMD 3 (`~/namd/namd3`, 10 CPU threads + GPU).
Two scripts in this folder do the work:

| Script | What it does |
|---|---|
| `run_queue.sh` | Runs **all** the simulations of the project one after another (see the list below). |
| `run_namd.sh` | Runs the steps of **one** `namd` folder (6.1 → 6.6, then 7). `run_queue.sh` uses it. |

## What the queue runs

| # | System | Folder | Temperature | Steps |
|---|---|---|---|---|
| 1 | 64 DMPC, 3200 water, 8 Na⁺ + 8 Cl⁻ | `charmm-N-64/charmm-gui-9562164596/namd` | 303.15 K | 6.1 → 7 |
| 2 | 128 DMPC, 6400 water, 17 Na⁺ + 17 Cl⁻ | `charmm-ORG/charmm-gui-9562151784/namd` | 303.15 K | 6.1 → 7 |
| 3 | 192 DMPC, 9600 water, 26 Na⁺ + 26 Cl⁻ | `charmm-N-192/charmm-gui-9562196952/namd` | 303.15 K | 6.1 → 7 |
| 4 | 128 DMPC (copy of #2) | `charmm-ORG/charmm-gui-9562151784/namd-TEMP-283` | 283 K | 6.6 → 7 |
| 5 | 128 DMPC (copy of #2) | `charmm-ORG/charmm-gui-9562151784/namd-TEMP-323` | 323 K | 6.6 → 7 |

Folders 4 and 5 are created automatically **after #2 has finished**. They are copies of
`charmm-ORG/.../namd` that keep its finished 6.1–6.5 (303.15 K) results; only `set temp`
in `step6.6_equilibration.inp` and `step7_production.inp` is changed to 283 or 323.

Expected time: roughly **15 hours in total** (13–17 h), based on the measured times of
earlier 128-lipid runs: about 2.3 h for N-64, 3.9 h for ORG, 5.8 h for N-192 and
1.6 h for each temperature.

## Before starting

1. **Plug in the laptop** and set it never to sleep while plugged in
   (PowerShell: `powercfg /change standby-timeout-ac 0`).
   If you will close the lid: *Control Panel → Power Options → Choose what closing the lid does →
   When plugged in: Do nothing*.
2. Put it on a hard, flat surface (ventilation). Earlier runs crashed the PC, most likely from heat or power.
3. Close heavy programs.

## Run everything

Open **Ubuntu** and run:

```bash
cd "/mnt/c/Users/Usuario/Documents/Arnau/Uni/4rt/SIMCOM/Github/Projecte Turri"
nohup bash run_queue.sh > run_queue.log 2>&1 &
tail -f run_queue.log
```

- `Ctrl+C` only stops *watching* the log; the simulations keep running.
- **Keep the Ubuntu window open.** If every WSL window is closed, WSL can shut down and kill NAMD.

## Check progress

```bash
cd "/mnt/c/Users/Usuario/Documents/Arnau/Uni/4rt/SIMCOM/Github/Projecte Turri"
tail -n 20 run_queue.log      # which system and step is running
pgrep -a namd3                # shows a line if NAMD is running, nothing if not
```

Each step writes a line like `Step 6.4 finished. WallClock: ...`.
The queue ends with `QUEUE FINISHED: all simulations done`.

## If it is interrupted (crash, PC closed, power cut)

Run **the same command again**:

```bash
cd "/mnt/c/Users/Usuario/Documents/Arnau/Uni/4rt/SIMCOM/Github/Projecte Turri"
nohup bash run_queue.sh > run_queue.log 2>&1 &
```

Steps that finished correctly are skipped (`Step X already finished, skipping`) and it
continues from the first unfinished step. The interrupted step starts again from its beginning.

If a step fails with an error, the queue stops there (`QUEUE STOPPED`), because the next
step needs its result. Look at the error in `run_queue.log` or in that step's log
(`stepX_equilibration.out`, or `output.log` for production).

## Stop everything

```bash
pkill -f run_queue.sh; pkill -f run_namd.sh; pkill namd3
```

## Run only one folder or some steps

```bash
cd "/mnt/c/Users/Usuario/Documents/Arnau/Uni/4rt/SIMCOM/Github/Projecte Turri"

bash run_namd.sh charmm-N-64/charmm-gui-9562164596/namd            # resume: runs what is not finished
bash run_namd.sh charmm-N-64/charmm-gui-9562164596/namd 6.6        # force: (re)run 6.6 and 7
bash run_namd.sh charmm-N-64/charmm-gui-9562164596/namd 6.4 6.6    # force: (re)run 6.4 to 6.6 only
```

Add `nohup ... > some.log 2>&1 &` around it to keep it running in the background.
Giving a first step **re-runs it even if it was finished**, overwriting its results.

## Where the results are

Each `namd` folder ends with:

| File | Content |
|---|---|
| `step6.X_equilibration.out` | log of each equilibration step |
| `output.log` | production log (energy, temperature, pressure → NAMDPlot) |
| `step7_production.dcd` | production trajectory, 200 frames (one every 2500 steps = 5 ps) → VMD |
| `step7_production.coor/.vel/.xsc` | final coordinates, velocities and box |

To load a trajectory in VMD: `step5_input.psf` + `step7_production.dcd` from the same folder.

## Analysis and plots (after the simulations)

From Windows (PowerShell or Git Bash, **not** WSL; it uses the Windows VMD install):

```bash
cd "C:/Users/Usuario/Documents/Arnau/Uni/4rt/SIMCOM/Github/Projecte Turri"
python Plots/run_analysis.py
```

For every finished run it writes energy, temperature, pressure vs pressure average,
g(r), N(r) and hydrogen bonds into `Plots/<Section>/` and rebuilds
`Plots/_output/Plots_report.pdf`. Runs that are not finished are skipped, so it can be
run again after each simulation ends. Sections: `Step3` (ORG, 303.15 K),
`Q1_T283K`, `Q1_T323K`, `Q2_N64`, `Q2_N192` (see `RUNS` in `Plots/run_analysis.py`).
VMD settings (selections, H-bond cutoffs) are in `Plots/vmd_analysis.tcl`.

## Technical notes

- NAMD runs in a copy of each folder on the WSL Linux disk (`~/namd-runs/...`), and the
  results are copied back to the Windows folder after every step. Writing directly to
  `/mnt/c` failed before with `Permission denied` while Windows was reading the files.
- A step only counts as finished if its log contains `WallClock:` and no `FATAL ERROR`.
  NAMD prints `End of program` even after a crash, so that line is not used.
- `step6.1_equilibration.inp` reads `../step5_assembly.str`, so a `namd` folder must stay
  inside its `charmm-gui-XXXXXXXXXX` folder.
- `dcdfreq` in `step7_production.inp` was changed from CHARMM-GUI's 50000 to **2500**
  (the same as in the earlier runs), so the production trajectories have 200 frames, not 10.
- To change the number of CPU threads, edit `NAMD_OPTS` in `run_namd.sh`.
- The old system (now archived in `v1/charmm-gui-DMPC-NaCl-membrane/`, see `v1/README.md`)
  contains **KCl** (`POT`), not NaCl. The three new systems contain NaCl (`SOD`).
