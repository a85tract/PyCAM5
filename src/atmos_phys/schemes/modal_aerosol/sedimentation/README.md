# Aerosol sedimentation process group

This directory contains the B06 CCPP-style boundaries for
`dust_sediment_tend_run` and `d3ddflux_run`.  The dust sedimentation body and
its private transport helpers move together; the original
`dust_sediment_mod` is a compatibility alias.  `d3ddflux` remains the CAM
adapter for the extracted array kernel because its host module also owns
other dry-deposition processes from B04 and B05.

Both processes have independent `ap_<process>_run` production timers.  Full
B06 acceptance remains gated on the integration lane's fresh 50-step BFB
comparison.
