continueflag = true

while (continueflag) {
	// Get background
	background = File.openDialog("Select Background");
	background_dir = File.getParent(background);
	background_filename = File.getName(background);
	// Get image stack
	image_stack = File.openDialog("Select image stack");
	image_stack_dir = File.getParent(image_stack);
	image_stack_filename = File.getName(image_stack);
	// Open them
	open(background);
	run("TIFF Virtual Stack...", "open=["+image_stack+"]");
	// Perform background subtraction
	imageCalculator("Subtract create 32-bit stack", image_stack_filename,background_filename);
	selectWindow("Result of " + image_stack_filename);
	run("8-bit");
	// Save background subtracted files in a different folder
	baseNameEnd=indexOf(image_stack_filename, ".tif"); 
	baseName=substring(image_stack_filename, 0, baseNameEnd);
	baseFolder = image_stack_dir + "/Processed/";
	File.makeDirectory(baseFolder); 	
	outfilename = baseFolder  + baseName + "_background_subtracted.tif";
	saveAs("Tiff", outfilename);
	
	// close all windows
	while (nImages>0) { 
	  selectImage(nImages); 
	  close(); 
	} 

	// Load the background subtracted image as a virtual stack
	run("TIFF Virtual Stack...", "open=["+outfilename+"]");
	// Obtain average intensity
	run("Z Project...", "projection=[Average Intensity]");
	outfilename = baseFolder  + baseName + "_background_subtracted_average.tif";
	saveAs("Tiff", outfilename);

    // close all windows
	while (nImages>0) { 
	  selectImage(nImages); 
	  close(); 
	} 	

	continueflag = getBoolean('Process one more?', 'Yes!', 'Nope');
}

