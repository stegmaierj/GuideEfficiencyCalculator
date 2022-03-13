# GuideEfficiencyCalculator

## Step 1: Upload Sequence Trace Files
		Click on the button named <i>Step 1: Upload Files</i> to switch to the file upload dialog. Drag&Drop the input files you want to process to the drop area or use the <i>Browse</i> button to open a file selection dialog.
		Valid input files should be in the <i>*.scf</i> format. Please convert input files in different formats, e.g., using software tools such as <a href="http://biologylabs.utah.edu/jorgensen/wayned/ape/" target="_blank">ApE</a> or using the command-line tool <a href="http://sourceforge.net/projects/staden/files/io_lib/1.12.2/io_lib-1.12.2.tar.gz" target="_blank">io lib</a>.
		Note that previous uploads will be overwritten, i.e., make sure to save your results before you upload another set of sequence files.
		
		<br><br><br>


## Step 2: Sequence Preview and Selection of the PAM Location
		After all files you want to process are uploaded successfully, you can proceed with the preview generation and the selection of the PAM location by clicking <i>Step 2: Preview / PAM Selection</i>. Note that this process may take a few seconds to minutes depending on the number of sequences you provided.
		This generates a preview image for each of the uploaded sequence files (see image below). Use the associated edit fields to enter the location of interest (<i>L</i>) for each of the sequences.
		The location of interest should be either set to the PAM location or to a location close to the expected cutting location.
		If the phred score of the read-out drops at a particular position, the location of the steepest quality drop can also be used. Be sure to specify this location separately for each of the sequences. 
		Clicking on one of the preview images will open a magnified view to be able to localize the location of interest more precisely.
		Finally, enter a positive integer value that specifies the analysis radius <i>R</i> (default value <i>R=50</i>). The guide efficiency calculations are based on comparing the <i>R</i> bases before and after the location of interest <i>L</i>, respectively.
		<br><br>
		<center>
		<img src="src/images/EfficiencyCalculationsOverview.png" width="100%">
		</center>
		
		<br><br><br>

## Step 3: Process Data and Download the Results
		Once all locations of interest are set properly and the analysis radius was specified you can continue with the actual efficiency prediction by clicking <i>Step 3: Process Data</i>. 
		This will automatically process all uploaded sequences with the specified parameters for <i>R</i> and <i>L</i> yielding overview plots of the calculations. Note that this process may take a few seconds to minutes depending on the number of sequences you provided.
		
		Use the download links on the right hand side of each of the preview images to save the results to disk. A zip archive containing all generated results can be downloaded at the top of the page. To proceed with another set of sequences, simply switch to <i>Step 1: Upload Files</i> and repeat steps 1-3.
		
		<br><br><br>
		<h3>Citation</h3>
		If you use our tool for your research, please cite the following publication: C. Etard, S. Joshi, J. Stegmaier, R. Mikut, U. Str&auml;hle: "Tracking of Indels by DEcomposition is a Simple and Effective Method to Assess Efficiency of Guide RNAs in Zebrafish", Zebrafish, 2017.
