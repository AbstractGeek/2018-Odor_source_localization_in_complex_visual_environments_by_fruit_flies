# @File(label="Directory of the images sequence", style="directory") images_sequence_dir
# @String(label="Image File Extension", required=false, value=".tif") image_extension
# @String(label="Output Filename", required=false, value=".") stack_name
# @String(label="Output Filename", required=false, value=".tif") output_name

import os
from ij import IJ
from pprint import pprint

from net.imagej.axis import Axes
from net.imglib2.view import Views

# Find image files
images_sequence_dir = str(images_sequence_dir)
print(images_sequence_dir)
fnames = []
for root, _, filenames in os.walk(images_sequence_dir):
    for fname in filenames:
        if fname.endswith(image_extension):
            fnames.append(os.path.join(root, fname))
fnames = sorted(fnames)
pprint(fnames)

if len(fnames) < 1:
    raise Exception("Not image files found in %s" % images_sequence_dir)

# Open and stack images
stack = []
for fname in fnames:
    image = IJ.openImage(fname)
    # then display it.
    image.show()

IJ.run("Images to Stack", "name="+str(stack_name)+" title=[] use")
IJ.saveAs("Tiff", os.path.join(images_sequence_dir,str(output_name)))
