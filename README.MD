# Files

## particleFilter.py

This file contains all of the particle filter code and associated plots.

## particle_cloud.py

This file contains all of the logic that handles the uncertainties for the
particle filter.

# ESU and ECEF

All of the final graphs are done in the ECEF space, which is a universal xyz
plane across the entire earth. However, to convert AR poses to ECEF we use an
intermediary plane called the ESU plane. In the eastupsouth json you will see
that there is associated data for it and you can use that for the conversions.

# TODOs

The major todo is change everything from staying in the ECEF space to being in
the ESU space. The ESU space would have an origin that you would determine.
