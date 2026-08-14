# Manifest format

Excel or CSV. Must have at least the following three (?) columns:
* tube_id
* sample_id
* project

Anything in the run that is not a control (e.g. a sample labeled "P|S|NC#") or not found in the manifest will be discarded. Other columns I don't care about - maybe add to the data container as metadata/obs/coldata