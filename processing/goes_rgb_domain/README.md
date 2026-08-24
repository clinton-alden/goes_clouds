# Generic GOES RGB domain workflow

`run_task.sh` downloads C02/C05/C13 ABI CONUS radiances for an inclusive date
slice, orthorectifies them using one locked/shared DEM, builds Zarr inputs and
daily RGB composites, and validates time coverage and geographic orientation.

The Mammoth configuration uses a 1 degree square centered on 37.63740 N,
119.03390 W: west/east -119.53390/-118.53390 and south/north
37.13740/38.13740. `prepare_mammoth_jobs.py` creates an 83-task manifest and
prints the qsub command. It does not submit unless `--submit` is supplied.
Outputs go into the established `/glade/derecho/scratch/cdalden/mammoth`
directory, using the distinct `mammoth_1deg` domain label. All tasks share
`static/SRTMGL3_mammoth_1deg.tif`; an inter-process lock guarantees that the
OpenTopography DEM request occurs only once even when tasks start together.
