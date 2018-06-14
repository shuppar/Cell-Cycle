%% Trying out clearing background for grey scale image from Image Analysts' generalised code. It's merely a copy of his code, with unwanted parts removed.


function correctedImage = BackgroundCorrect(inputImage, nonuniformBackgroundImage)
try
	
	correctedImage = inputImage;
    		 modeledBackgroundImage = zeros(size(nonuniformBackgroundImage));
   		 noiselessImage = nonuniformBackgroundImage;
			maxValue = max(max(noiselessImage));
			modeledBackgroundImage = noiselessImage / maxValue;
			correctedImage = (inputImage ./ modeledBackgroundImage);
        
end
return;
