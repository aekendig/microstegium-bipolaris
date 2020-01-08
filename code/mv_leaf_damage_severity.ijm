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

input=//enter input folder path;
outputIm=//enter output image folder path;
outputRes=//enter output text results folder path;

suffix1=".tif"; //Store potential suffixes as variables
suffix2=".tiff";

processFolder(input);

// saves results for all images in a single file (moved from function processFolder) 
selectWindow("Summary"); 
//saveAs("Results", outputRes + "/All_Results_Mv_tif.tsv"); //change for appropriate suffix 
//saveAs("Results", outputRes + "/All_Results_Mv_tiff.tsv");
saveAs("Results", outputRes + "/All_Results_Mv_tiff_edited.tsv");

// function to scan folders/subfolders/files to find files with either correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix1))
			processFile(input, outputIm, outputRes, list[i]);
		if(endsWith(list[i], suffix2))
			processFile(input, outputIm, outputRes, list[i]);	
	}
 
}

// full process function
function processFile(input, outputIm, outputRes, file) {
	setBatchMode(true); // prevent image windows from opening while the script is running
	print("Processing: " + input + File.separator + file); // updates on progress
	open(input + "/" + file); // open file
	run("Set Scale...", "distance=0 global"); // measure in pixes
  	id = getTitle(); // get original image id
	RGBpic=getTitleStripExtension(); // get image name without extension
	run("Duplicate...", "title=leaf:"+RGBpic); // create a copy to segment leaf from background
	// color threshold function (HSB) to select the leaf and ignore the background
	min=newArray(3);
	max=newArray(3);
	filter=newArray(3);
	a=getTitle();
	run("HSB Stack");
	run("Convert Stack to Images");
	selectWindow("Hue");
	rename("0");
	selectWindow("Saturation");
	rename("1");
	selectWindow("Brightness");
	rename("2");
	min[0]=0;
	max[0]=73; // hue is 0 to 73
	filter[0]="pass";
	min[1]=3;
	max[1]=255; // saturation is 3 to 255
	filter[1]="pass";
	min[2]=20;
	max[2]=245; // brightness is 0 to 245
	filter[2]="pass";
	for (i=0;i<3;i++){
	  selectWindow(""+i);
	  setThreshold(min[i], max[i]);
	  run("Convert to Mask");
	  if (filter[i]=="stop")  run("Invert");
	}
	imageCalculator("AND create", "0","1");
	imageCalculator("AND create", "Result of 0","2");
	for (i=0;i<3;i++){
	  selectWindow(""+i);
	  close();
	}
	selectWindow("Result of 0");
	close();
	selectWindow("Result of Result of 0");
	rename(a);
	// end of color threshold function
	run("Fill Holes"); // edit ROI to include the leaf, but remove shadows on edges and dirt specks
	run("Close-");
	run("Fill Holes");
	run("Remove Outliers...", "radius=10 threshold=50 which=Bright");
	run("Analyze Particles...", "size=1000-Infinity pixel display include summarize add"); // measure leaf area
	selectWindow(id); // select original image
	roiManager("Show All"); // overlay the leaf ROIs made above
	if (roiManager("count") > 1) roiManager("Combine"); // combine all ROIs into one object if more than one
	if(roiManager("count") == 1) roiManager("Select", 0); // select one if only one ROI
	run("Clear Outside"); // clear the outside (blue background)
	run("Hide Overlay"); // remove ROI from image
	roiManager("Delete"); // clear ROI Manager for next image
	run("Duplicate...", "title=lesions:"+RGBpic); // picture to isolate lesions from leaf
	run("Duplicate...", "title=output:"+RGBpic); // picture to save lesion image
	run("Duplicate...", "title=greens:"+RGBpic); // picture to isolate green space from lesions
	selectWindow("lesions:"+RGBpic);
	// color threshold function (Lab) to select lesions and ignore leaf
	min=newArray(3);
	max=newArray(3);
	filter=newArray(3);
	a=getTitle();
	call("ij.plugin.frame.ColorThresholder.RGBtoLab"); // note RGB to Lab
	run("RGB Stack");
	run("Convert Stack to Images");
	selectWindow("Red");
	rename("0");
	selectWindow("Green");
	rename("1");
	selectWindow("Blue");
	rename("2");
	min[0]=0;
	max[0]=253; // L* is 0 to 253
	filter[0]="pass";
	min[1]=125;
	max[1]=169; // a* is 125 to 169
	filter[1]="pass";
	min[2]=0;
	max[2]=255; // b* is 0 to 255
	filter[2]="pass";
	for (i=0;i<3;i++){
	  selectWindow(""+i);
	  setThreshold(min[i], max[i]);
	  run("Convert to Mask");
	  if (filter[i]=="stop")  run("Invert");
	}
	imageCalculator("AND create", "0","1");
	imageCalculator("AND create", "Result of 0","2");
	for (i=0;i<3;i++){
	  selectWindow(""+i);
	  close();
	}
	selectWindow("Result of 0");
	close();
	selectWindow("Result of Result of 0");
	rename(a);
	// end of color threshold function
	run("Close-"); // edit ROI to include full lesions
	run("Fill Holes");
	run("Dilate");
	run("Dilate");
	run("Close-");
	run("Fill Holes");
	run("Erode");
	run("Erode");
	run("Remove Outliers...", "radius=5 threshold=50 which=Dark");
	run("Analyze Particles...", "size=0-Infinity pixel display include summarize add"); // measure lesion area
	run("Colors...", "foreground=white background=white selection=yellow"); // set new line color to white
	selectWindow("output:"+RGBpic); // select image of leaf
	roiManager("Show All"); // overlay the lesion ROIs made above
	roiManager("Draw"); // draw lesions on it
	selectWindow("greens:"+RGBpic); // select the image to separate the lesions from the green space
	roiManager("Show All"); // overlay the lesion ROIs made above

	//Finish the process and save the results if no lesions
	if(roiManager("count") < 1){
		selectWindow("output:"+RGBpic); // select image of leaf
		roiManager("Show All"); // overlay the green ROIs made above
		roiManager("Draw"); // draw green area on it
		save(outputIm + "/" + RGBpic + "_lesions" + ".tif"); // save image with lesions and green area
		if (roiManager("count") > 0) roiManager("Delete"); // clear ROI Manager for next image
		selectWindow("Results"); 
		saveAs("Results", outputRes + "/" + RGBpic + "_Results.tsv"); // save info on lesions and green area
		run("Clear Results");
		run("Close All");
		run("Collect Garbage");
	}

	//Remove green space in lesions before saving results if lesions are present
	if(roiManager("count") > 0){
		if (roiManager("count") > 1) roiManager("Combine"); // combine all ROIs into one object if more than one
		if(roiManager("count") == 1) roiManager("Select", 0); // select one if only one ROI
		run("Clear Outside"); // clear leaf outside lesions
		run("Hide Overlay"); // remove ROIs
		roiManager("Delete"); // clear ROI Manager for next image
		// color threshold function (HSB) to select the greens and ignore the lesions
		min=newArray(3);
		max=newArray(3);
		filter=newArray(3);
		a=getTitle();
		run("HSB Stack");
		run("Convert Stack to Images");
		selectWindow("Hue");
		rename("0");
		selectWindow("Saturation");
		rename("1");
		selectWindow("Brightness");
		rename("2");
		min[0]=38;
		max[0]=120; // hue is 38 to 120
		filter[0]="pass";
		min[1]=75;
		max[1]=255; // saturation is 75 to 255
		filter[1]="pass";
		min[2]=0;
		max[2]=254; // brightness is 0 to 254
		filter[2]="pass";
		for (i=0;i<3;i++){
		  selectWindow(""+i);
		  setThreshold(min[i], max[i]);
		  run("Convert to Mask");
		  if (filter[i]=="stop")  run("Invert");
		}
		imageCalculator("AND create", "0","1");
		imageCalculator("AND create", "Result of 0","2");
		for (i=0;i<3;i++){
		  selectWindow(""+i);
		  close();
		}
		selectWindow("Result of 0");
		close();
		selectWindow("Result of Result of 0");
		rename(a);
		// end of color threshold function
		run("Despeckle"); // refine image
		run("Dilate");
		run("Dilate");
		run("Close-");
		run("Erode");
		run("Erode");
		run("Remove Outliers...", "radius=5 threshold=50 which=Dark");
		run("Analyze Particles...", "display include summarize add"); // measure area of greens
		run("Colors...", "foreground=cyan background=white selection=yellow"); // set new line color to blue
		selectWindow("output:"+RGBpic); // select image of leaf
		roiManager("Show All"); // overlay the green ROIs made above
		roiManager("Draw"); // draw green area on it
		save(outputIm + "/" + RGBpic + "_lesions" + ".tif"); // save image with lesions and green area
		if (roiManager("count") > 0) roiManager("Delete"); // clear ROI Manager for next image
		selectWindow("Results"); 
		saveAs("Results", outputRes + "/" + RGBpic + "_Results.tsv"); // save info on lesions and green area
		run("Clear Results");
		run("Close All");
		run("Collect Garbage");
	}
	
}



