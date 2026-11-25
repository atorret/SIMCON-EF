# vmd_analysis.tcl - structure analysis of one NAMD production run (VMD, text mode).
#
# Not run directly: run_analysis.py writes a small header that sets
#   PSF   structure file (step5_input.psf)
#   DCD   trajectory (step7_production.dcd)
#   OUT   output folder (Plots/<Section>)
#   STEP  frame stride (1 = every frame)
#   LAST  last trajectory frame to read (-1 = all; set lower when the end is corrupt)
# and then sources this file with  vmd -dispdev text -e <header>.
#
# Writes into OUT:
#   G(r)/<pair>.dat               r  g(r)
#   N(r)/<pair>.dat               r  n(r)   (running coordination number)
#   Hydrogen_bonds.dat            frame  number of water-water H-bonds
#   Hydrogen_bonds_per_water.dat  frame  H-bonds per water molecule (2 N_HB / N_water)
#
# Selections and settings follow the project brief.

# g(r) pairs: file name, selection 1, selection 2
set PAIRS {
    Ow_Hw  {type OT}   {type HT}
    Ow_Na  {type OT}   {type SOD}
    Hw_Cl  {type HT}   {type CLA}
    N_Na   {type NTL}  {type SOD}
    N_Cl   {type NTL}  {type CLA}
    P_Na   {type PL}   {type SOD}
    P_Cl   {type PL}   {type CLA}
}
set GOFR_DELTA 0.1     ;# bin width (Å)
set GOFR_RMAX  10.0    ;# maximum r (Å)
set HB_SEL     {resname TIP3}
set HB_DIST    3.6     ;# donor-acceptor distance (Å)
set HB_ANGLE   30      ;# angle cutoff (degrees)

proc say {msg} { puts "ANALYSIS) $msg" }

proc write_columns {path xs ys} {
    set fh [open $path w]
    foreach x $xs y $ys { puts $fh "$x $y" }
    close $fh
}

if {![info exists STEP]} { set STEP 1 }
if {![info exists LAST]} { set LAST -1 }

mol new $PSF type psf waitfor all
mol addfile $DCD type dcd first 0 last $LAST step $STEP waitfor all
set nframes [molinfo top get numframes]
say "loaded $nframes frames (stride $STEP)"
if {$nframes < 1} { error "no frames in $DCD" }

file mkdir [file join $OUT "G(r)"] [file join $OUT "N(r)"]

# ---- Radial distribution functions and coordination numbers ----------------
foreach {name s1 s2} $PAIRS {
    set a [atomselect top $s1]
    set b [atomselect top $s2]
    set na [$a num]
    set nb [$b num]
    if {$na == 0 || $nb == 0} {
        say "WARNING $name skipped: empty selection ('$s1' = $na atoms, '$s2' = $nb atoms)"
        $a delete; $b delete
        continue
    }
    set t0 [clock milliseconds]
    # result: {r} {g(r)} {number integral} {histogram} {frames total skipped processed}
    set res [measure gofr $a $b delta $GOFR_DELTA rmax $GOFR_RMAX usepbc 1 \
                 first 0 last -1 step 1]
    lassign $res r g n
    write_columns [file join $OUT "G(r)" "$name.dat"] $r $g
    write_columns [file join $OUT "N(r)" "$name.dat"] $r $n
    say [format "%-6s %5d x %5d atoms  %.1f s" $name $na $nb \
             [expr {([clock milliseconds] - $t0) / 1000.0}]]
    $a delete; $b delete
}

# ---- Hydrogen bonds -----------------------------------------------------------
package require hbonds
set water [atomselect top $HB_SEL]
set nwater [llength [lsort -unique [$water get residue]]]
set hbfile [file join $OUT Hydrogen_bonds.dat]
set t0 [clock milliseconds]
hbonds -sel1 $water -dist $HB_DIST -ang $HB_ANGLE -frames all -plot no \
       -writefile yes -outfile $hbfile -log [file join $OUT hbonds_vmd.log]
file delete [file join $OUT hbonds_vmd.log]
say [format "H-bonds  %d water molecules  %.1f s" $nwater \
         [expr {([clock milliseconds] - $t0) / 1000.0}]]

# Each H-bond joins two molecules, so per molecule it counts twice
# (comparable to ~4 per molecule in bulk water).
set fin [open $hbfile r]
set fout [open [file join $OUT Hydrogen_bonds_per_water.dat] w]
while {[gets $fin line] >= 0} {
    if {[llength $line] < 2} continue
    lassign $line frame count
    puts $fout [format "%d %.4f" $frame [expr {2.0 * $count / $nwater}]]
}
close $fin
close $fout
$water delete

say "DONE"
# (the header written by run_analysis.py calls quit)
