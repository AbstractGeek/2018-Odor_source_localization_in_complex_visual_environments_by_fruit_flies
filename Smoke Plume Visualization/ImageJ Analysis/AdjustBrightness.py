# @File(label="Directory of the images sequence", style="directory") images_sequence_dir
# @String(label="Image File Extension", required=false, value=".tif") image_extension
# @File(label="Directory for the output images", style="directory") output_dir

import os
from ij import IJ
from pprint import pprint

from net.imagej.axis import Axes
from net.imglib2.view import Views

# Find image files
images_sequence_dir = str(images_sequence_dir)
print(images_sequence_dir)
fnames = []
filenames = []
for fname in os.listdir(images_sequence_dir):    
        if fname.endswith(image_extension):
            fnames.append(os.path.join(images_sequence_dir, fname))
            filenames.append(os.path.join(str(output_dir),fname))
fnames = sorted(fnames)
pprint(fnames)
pprint(filenames)

if len(fnames) < 1:
    raise Exception("Not image files found in %s" % images_sequence_dir)

# Add the main function below
for i, fname in enumerate(fnames):
    image = IJ.openImage(fname)
    # then display it.
    image.show()
    IJ.run("Enhance Contrast")
    IJ.run("Apply LUT")
    IJ.saveAs("PNG", os.path.join(filenames[i]))
    IJ.run("Close")
