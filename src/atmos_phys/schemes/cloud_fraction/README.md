# Cloud-fraction process group

This directory contains the B06 `cloud_fraction_fice_run` numerical process
and its CCPP metadata.  The original `cldfrc_fice` routine remains as a thin
CAM adapter, preserving its public interface and owning the
`ap_cloud_fraction_fice_run` production timer.

The later B07 cloud-fraction extraction extends this group; this B06 boundary
is limited to the temperature-based liquid/ice/snow partition.
