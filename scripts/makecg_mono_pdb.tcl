#script to create a trajectory of CG beads for desired atom selections
mol new GGGGG.pdb type pdf 

set outfile GGGGG.cg.pdb
set rgyroutfile GGGGG.cg.rgyr

#enter the atomselections in selList
set selList {"resid 1 and name N" "resid 1 and name H" "resid 1 and name CA" "resid 1 and name HA2" "resid 1 and name HA3" "resid 1 and name C" "resid 1 and name O" "resid 1 and name OXT" "resid 2 and name N" "resid 2 and name H" "resid 2 and name CA" "resid 2 and name HA2" "resid 2 and name HA3" "resid 2 and name C" "resid 2 and name O" "resid 3 and name N" "resid 3 and name H" "resid 3 and name CA" "resid 3 and name HA2" "resid 3 and name HA3" "resid 3 and name C" "resid 3 and name O" "resid 4 and name N" "resid 4 and name H" "resid 4 and name CA" "resid 4 and name HA2" "resid 4 and name HA3" "resid 4 and name C" "resid 4 and name O" "resid 5 and name NA" "resid 5 and name C8" "resid 5 and name O1" "resid 5 and name C9" "resid 5 and name O2" "resid 5 and name C10" "resid 5 and name C11" "resid 5 and name H16" "resid 5 and name C12" "resid 5 and name C13" "resid 5 and name C14" "resid 5 and name C15" "resid 5 and name C16" "resid 5 and name H17" "resid 5 and name C17" "resid 5 and name H18" "resid 5 and name C18" "resid 5 and name C19" "resid 5 and name C20" "resid 5 and name C21" "resid 5 and name H19" "resid 5 and name H20" "resid 5 and name C22" "resid 5 and name C23" "resid 5 and name C24" "resid 5 and name H21" "resid 5 and name C25" "resid 5 and name H22" "resid 5 and name C26" "resid 5 and name C27" "resid 5 and name O3" "resid 5 and name C28" "resid 5 and name C29" "resid 5 and name C30" "resid 5 and name O4" "resid 5 and name NB" "resid 5 and name C31" "resid 5 and name H23" "resid 5 and name CA" "resid 5 and name OA" "resid 5 and name CA1" "resid 5 and name HA1" "resid 5 and name HA2" "resid 5 and name CB1" "resid 5 and name HB1" "resid 5 and name HB2" "resid 5 and name CB" "resid 5 and name OB" "resid 6 and name N" "resid 6 and name H" "resid 6 and name CA" "resid 6 and name HA2" "resid 6 and name HA3" "resid 6 and name C" "resid 6 and name O" "resid 7 and name N" "resid 7 and name H" "resid 7 and name CA" "resid 7 and name HA2" "resid 7 and name HA3" "resid 7 and name C" "resid 7 and name O" "resid 8 and name N" "resid 8 and name H" "resid 8 and name CA" "resid 8 and name HA2" "resid 8 and name HA3" "resid 8 and name C" "resid 8 and name O" "resid 9 and name N" "resid 9 and name H" "resid 9 and name CA" "resid 9 and name HA2" "resid 9 and name HA3" "resid 9 and name C" "resid 9 and name O" "resid 9 and name OXT"}

#enter the atom names of the CG particles
set CGNAMElist { "N" "H" "CA" "HA2" "HA3" "C" "O" "OXT" "N" "H" "CA" "HA2" "HA3" "C" "O" "N" "H" "CA" "HA2" "HA3" "C" "O" "N" "H" "CA" "HA2" "HA3" "C" "O" "NA" "C8" "O1" "C9" "O2" "C10" "C11" "H16" "C12" "C13" "C14" "C15" "C16" "H17" "C17" "H18" "C18" "C19" "C20" "C21" "H19" "H20" "C22" "C23" "C24" "H21" "C25" "H22" "C26" "C27" "O3" "C28" "C29" "C30" "O4" "NB" "C31" "H23" "CAT" "OAT" "CA1" "HA1" "HA2" "CB1" "HB1" "HB2" "CBT" "OBT" "N" "H" "CA" "HA2" "HA3" "C" "O" "N" "H" "CA" "HA2" "HA3" "C" "O" "N" "H" "CA" "HA2" "HA3" "C" "O" "N" "H" "CA" "HA2" "HA3" "C" "O" "OXT" } 

set CGRESNAMElist { "GC1" "GC1" "GC1" "GC1" "GC1" "GC1" "GC1" "GC1" "G1" "G1" "G1" "G1" "G1" "G1" "G1" "G1" "G1" "G1" "G1" "G1" "G1" "G1" "G1" "G1" "G1" "G1" "G1" "G1" "G1" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI" "PDI"  "G2" "G2" "G2" "G2" "G2" "G2" "G2" "G2" "G2" "G2" "G2" "G2" "G2" "G2" "G2" "G2" "G2" "G2" "G2" "G2" "G2" "GC2" "GC2" "GC2" "GC2" "GC2" "GC2" "GC2" "GC2"}

set RESIDlist {"1" "1" "1" "1" "1" "1" "1" "1" "2" "2" "2" "2" "2" "2" "2" "3" "3" "3" "3" "3" "3" "3" "4" "4" "4" "4" "4" "4" "4" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "5" "6" "6" "6" "6" "6" "6" "6" "7" "7" "7" "7" "7" "7" "7" "8" "8" "8" "8" "8" "8" "8" "9" "9" "9" "9" "9" "9" "9" "9"}

set RESNAME "RNA"

set CHAINNAME "A"

#create a pdb file name for the CG particle output
set outstream [open $outfile w]

#fine the number of selections and the number of frames in the loaded trajectory
set selListLength [llength $selList]
set nFrames [molinfo top get numframes]

#create atom selections for each of the selection above in selList by creating an array of atom selections
for {set j 0} {$j < $selListLength} {incr j} {
	set sel($j) [atomselect top [lindex $selList $j]]
	set srgyr($j) [expr 0]
}


#for each frame, update the selection to that frame, find the x y and z positions of the COM, and output the information into the PDB format
for {set k 0} {$k < $nFrames} {incr k} {
	for {set j 0} {$j < $selListLength} {incr j} {
		$sel($j) frame $k
		$sel($j) update
		set scom [measure center $sel($j) weight mass]
		set scomx [lindex $scom 0]
		set scomy [lindex $scom 1]
		set scomz [lindex $scom 2]
	    puts $outstream [format "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" $j [lindex $CGNAMElist $j] [lindex $CGRESNAMElist $j] $CHAINNAME [lindex $RESIDlist $j] $scomx $scomy $scomz 1.00 0.00]
		flush $outstream
		set srgyr($j) [expr $srgyr($j)+[measure rgyr $sel($j)]]
	}
	puts "Frame $k finished. Next frame."
	puts $outstream "END"
	flush $outstream
}

#close file
close $outstream


set rgyrout [open $rgyroutfile w]
for {set j 0} {$j < $selListLength} {incr j} {
	set srgyr($j) [expr $srgyr($j)/$nFrames]
	puts $rgyrout [format "%4d %8.3f" [expr $j]  $srgyr($j)]
	flush $rgyrout
}

close $rgyrout
	
## make a dcd out of pdb
#mol new rna.50bps.v01.cg.psf type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
#mol addfile $outfile type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
#set num [molinfo top get numframes] 
#incr num -1
#animate write dcd rna.50bps.v01.cg.dcd beg 0 end $num skip 1 waitfor all 
#animate delete all
#quit
