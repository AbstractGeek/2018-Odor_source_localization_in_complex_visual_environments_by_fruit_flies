continueflag = true

while (continueflag) {
	// Get image stack
	image_stack = File.openDialog("Select image stack");
	image_stack_dir = File.getParent(image_stack);
	image_stack_filename = File.getName(image_stack);
	// Open them
	run("TIFF Virtual Stack...", "open=["+image_stack+"]");
	// Obtain average intensity
	run("Z Project...", "projection=[Average Intensity]");
	
	// Save background subtracted files in a different folder
	baseNameEnd=indexOf(image_stack_filename, ".tif"); 
	baseName=substring(image_stack_filename, 0, baseNameEnd);
	baseFolder = image_stack_dir + "/Processed/";
	File.makeDirectory(baseFolder); 	
	outfilename = baseFolder  + baseName + "_average.tif";
	saveAs("Tiff", outfilename);
	
	// close all windows
	while (nImages>0) { 
	  selectImage(nImages); 
	  close(); 
	} 

	continueflag = getBoolean('Process one more?', 'Yes!', 'Nope');
}

