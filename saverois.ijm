count = roiManager("count");
for (i = 0; i < count; i++) {
	roiManager("select", i);
	Roi.getCoordinates(xx, yy);
	Array.print(xx);
	Array.print(yy);	
}

