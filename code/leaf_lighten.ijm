// go to Options...
// make sure "Black Background" is not selected


// strip extension off image title
// will leaf "f" if it takes off .tif first and it's a .tiff
function getTitleStripExtension() { 
  t = getTitle(); 
  t = replace(t, ".tiff", "");   
  t = replace(t, ".tif", "");             
  t = replace(t, ".lif", "");       
  t = replace(t, ".lsm", "");     
  t = replace(t, ".czi", "");       
  t = replace(t, ".nd2", "");     
  return t; 
} 

//enter home folder (errors in path names will result in the error "No window with the title 'Summary' found"
home="/Users/AmyKendig/Dropbox (UFL)/big-oaks-field-experiment-2018-2019/leaf-scans/leaf-scans-add-on";
//enter input folder path
input=home + "/scans/mv"; 
//enter output image folder path
outputIm=home + "/lightened-scans";

suffix1=".tiff"; //Store potential suffixes as variables
suffix2=".tif";

processFolder(input);

// function to scan folders/subfolders/files to find files with either correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix1))
			processFile(input, outputIm, list[i]);
		if(endsWith(list[i], suffix2))
			processFile(input, outputIm, list[i]);	
	}
}

// full process function
function processFile(input, outputIm, file) {
	setBatchMode(true); // prevent image windows from opening while the script is running
	print("Processing: " + input + File.separator + file); // updates on progress
	open(input + "/" + file); // open file
  	id = getTitle(); // get original image id
	RGBpic=getTitleStripExtension(); // get image name without extension
	run("Duplicate...", "title=light:"+RGBpic); // create a copy to segment leaf from background
	//run("Brightness/Contrast...");
	resetMinAndMax();
	setMinAndMax(-36, 218);
	run("Apply LUT");
	//save the image
	selectWindow("light:"+RGBpic); // select image of leaf
	save(outputIm + "/" + RGBpic + ".tif"); // save image with original name
	run("Close All");
	run("Collect Garbage");
}