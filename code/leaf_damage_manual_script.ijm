// go to Options...
// make sure "Black Background" is not selected

suffix1=".tif"; //Store potential suffixes as variables
suffix2=".tiff";

	run("Set Scale...", "distance=0 global"); // measure in pixels
  	id = getTitle(); // get original image id
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
	min[1]=127;
	max[1]=169; // a* is 127 to 169
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
	// run("Close-"); // seems to connect a lot of the leaf background
	run("Remove Outliers...", "radius=5 threshold=50 which=Dark"); // added this in for script 3
	run("Close-");
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

	//Save the results if no lesions
	save(outputIm + "/" + RGBpic + "_lesions" + ".tif"); // save image with lesions and green area
	if (roiManager("count") > 0) roiManager("Delete"); // clear ROI Manager for next image
	selectWindow("Results"); 
	saveAs("Results", outputRes + "/" + RGBpic + "_Results.tsv"); // save info on lesions and green area
	run("Clear Results");
	run("Close All");
	run("Collect Garbage");