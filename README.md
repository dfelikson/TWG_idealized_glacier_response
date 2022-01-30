# TWG_idealized_glacier_response

This repository contains scripts to configure the Ice Sheet System Model to perform the simulations described in the manuscript *Seasonal tidewater glacier terminus oscillations bias multi-decadal projections of ice mass change*, submitted for review to the Journal of Geophysical Research in May 2021. The repository contains the following files:

1. **runme.m**: Main script that creates mesh, parameterizes the model, and runs all transient simulations. Different parts of this script can be run by specifying which step(s) are executed in the first line.
1. **move_terminus_levelset.m**: Moves the zero-contour of a levelset by a specified distance. This allows the user to specify retreat or advance of the glacier terminus by some distance.
1. **MISMIP_TWG.par**: Script containing model parameters. Step 2 in the runme.m script uses this file to establish model parameters.
1. **Exp/Domain_TWG.exp**: This file contains the coordinates that specify the outline of the model domain. The file is in ARGUS format and can be read with an ASCII text reader.
1. **Exp/Front_TWG.exp**: This file contains the coordinates that outline the model mesh vertices that are located at the ice front (i.e., the ice-ocean boundary). The file is in ARGUS format and can be read with an ASCII text reader.

The runme.m script saves model output to the **Models_TWG/** directory.
