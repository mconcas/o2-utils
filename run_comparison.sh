#!/bin/bash

# run reference
root -l -q run_trac_ca_its.C++\(false\) && mv dbg_ITSTrackerCPU.root nosmooth/dbg_ITSTrackerCPU.root

# run smoother
root -l -q run_trac_ca_its.C++\(true\) && mv dbg_ITSTrackerCPU.root smooth/dbg_ITSTrackerCPU.root

# run comparison
root -l compareTrackerResults.C++