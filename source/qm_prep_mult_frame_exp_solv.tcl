proc geom_center {selection} {
        # set the geometrical center to 0
	#set list_center {}
        set gc [veczero]
        # [$selection get {x y z}] returns a list of {x y z} 
        #    values (one per atoms) so get each term one by one
        foreach coord [$selection get {x y z}] {
           # sum up the coordinates
           set gc [vecadd $gc $coord]
        }
        # and scale by the inverse of the number of atoms
	#set center [vecscale [expr 1.0 /[$selection num]] $gc]
	return [vecscale [expr 1.0 /[$selection num]] $gc]
	#lappend list_center $center
        #return $list_center
}

proc findmax { items } {
	set max [lindex $items 0]
        foreach i $items { 
		if { $i > $max } {
			set max $i
                }
        }
	return $max
}

proc sortDictByValue {dict args} {
    set lst {}
    dict for {k v} $dict {lappend lst [list $k $v]}
    return [concat {*}[lsort -real -index 1 {*}$args $lst]]
}

proc ldiff {list1 list2} {
	set result {}
	foreach item $list1 {if {[lsearch $list2 $item]==-1} {lappend result $item}}
	set result
}


set list_1 {}
set list_2 {}
for {set i 0} {$i < 104} {incr i} {
        set other_filename I1_all-atom_fixed_$i-ns.pdb
        mol load pdb $other_filename
	set other_sel_water [atomselect top "water and within 10 of (resname CU1P DTT DTBD)"]
	set water_all_resid [$other_sel_water get resid]
	#puts $water_all_resid
	set water_unique_resid [lsort -unique $water_all_resid]
	#puts $water_unique_resid
	set len_water [llength $water_unique_resid]
	append list_1 $len_water ","
	append list_2 "\t"  $len_water
	mol delete $i
}
set other_filename I1_all-atom_fixed_$i-ns.pdb
mol load pdb $other_filename
set other_sel_water [atomselect top "water and within 10 of (resname CU1P DTT DTBD)"]
set water_all_resid [$other_sel_water get resid]
#puts $water_all_resid
set water_unique_resid [lsort -unique $water_all_resid]
#puts $water_unique_resid
set len_water [llength $water_unique_resid]
append list_1 $len_water
append list_2 "\t" $len_water
mol delete [expr $i]

#puts $list_1
#set base_water_number [expr min($list_1) ]
set base_water_number 50
#puts "water $base_water_number"
set base_water_index [lsearch -all -exact -integer $list_2 $base_water_number]
puts "index $base_water_index"
#puts [lindex $list_2  4]

set list_0 {}
for {set i 0} {$i < 104} {incr i} {
	lappend list_0 $i
}

set mod_list_0 [ldiff $list_0 $base_water_index]

foreach i $mod_list_0 {
		set other_filename I1_all-atom_fixed_$i-ns.pdb
        	mol load pdb $other_filename
        	set other_sel_water [atomselect top "water and within 10 of (resname CU1P DTT DTBD)"]
        	set water_all_resid [$other_sel_water get resid]
		set pol_sel_water [atomselect top "water and within 3.5 of name S1 S2 S3 S4 N2 O1 O2 CU1P"]
		set pol_sel_water_resid [$pol_sel_water get resid]
		set unique_pol_sel_water_resid [lsort -unique $pol_sel_water_resid]
		puts [llength $unique_pol_sel_water_resid]
        	set water_unique_resid [lsort -unique $water_all_resid]
		puts [llength $water_unique_resid]
		set mod_water_unique_resid [ldiff $water_unique_resid $unique_pol_sel_water_resid]
		puts [llength $mod_water_unique_resid]
        	#puts $water_unique_resid
		set list_1 {}
		set solute [atomselect top "resname CU1P DTT DTBD"]
        	set solute_center [geom_center $solute]
		puts $solute_center
		foreach unique_resid $mod_water_unique_resid {
			set solv [atomselect top "resid $unique_resid and name OH2"]
			set solv_coord [$solv get {x y z}]
			puts $solv_coord
			set dist [veclength [vecsub $solute_center [expr ($solv_coord)]]]
			lappend list_1 $dist
			dict set distances $unique_resid $dist
			puts $distances
		}

		set total_water [lindex $list_2 $i]
		puts $total_water
		set del_water [expr $total_water - $base_water_number]
		puts $del_water
		#puts $distances
		set distances_sorted [sortDictByValue $distances]
		#puts $distances_sorted
		#puts [dict size $distances_sorted]
		set new_distances [dict keys $distances_sorted]
		set mod_wat_resid [lrange $new_distances 0 end-$del_water]
		puts [llength $mod_wat_resid]
		set final_wat_resid [concat $unique_pol_sel_water_resid $mod_wat_resid]
		set water_no [llength $final_wat_resid]
		puts "water $water_no"
		set final_selection [atomselect top "(resname CU1P DTT DTBD) or resid $final_wat_resid"]
		$final_selection writepdb qm_system_I1_$i-ns.pdb
		set distances {}
		set distances_sorted {}
		#puts $distances_$i
}
foreach i $base_water_index {
		set other_filename I1_all-atom_fixed_$i-ns.pdb
                mol load pdb $other_filename
		set other_sel_water [atomselect top "water and within 10 of (resname CU1P DTT DTBD)"]
                set water_all_resid [$other_sel_water get resid]
		set water_unique_resid [lsort -unique $water_all_resid]
		set final_selection [atomselect top "(resname CU1P DTT DTBD) or resid $water_unique_resid"]
		$final_selection writepdb qm_system_I1_$i-ns.pdb

}
exit	
