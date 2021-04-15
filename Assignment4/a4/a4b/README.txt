CISC 472: Assignment 3
Evelyn Yach (20071956)
Daniel Oh (20063998)

5. For the CT/MR measures, comment very briefly (two sentences for each, at most) on each of:
     - the effect of rotation and translation away from the perfect registration
     - the smoothness of the MI, RMS, and NCC measures
     - the capture range of MI, RMS, and NCC

The effect of rotation and translation away from perfect registration seems to smooth out the plots. Sharp maxima 
and minima, seem to even out becoming less extreme.

MI measure seems to be the smoothest, followed by NCC and RMS which seem to have relatively the same smoothness.

Examining the axes, NCC seems to have the largest capture range, followed by RMS and then MI which is the smallest.


6. For CT/CT and CT/MR, why do RMS and NCC get larger with very large translations?  Why does MI *not* get larger?

We see RMS and NCC get larger with large translations and not MI because when there is only a small overlap between 
the two images RMS and NCC continue to make calculations in each of their respective methods. However, since MI relies 
on a predictive method between each pixel on one image and each on the other, the negative space (which gives limited 
info for predictions) causes a decline in the plot. This phenomena is illustrated on the graphs as RMS and NCC being 
more bowl shaped towards the edges, while MI has more of a dome appearance. 


7. What is your recommendation for a similarity measure for registering CT to MR images?  Briefly justify your answer.

Considering that both images will be coming from different imaging modalities, this means that it will be unlikely 
that the pixels of the two images will not have the same normalized values at optimal alignment. MI should be used 
as it does not require the images to be from the same modalities or using the same parameters. The intensity for 
one pixel in an image does not directly correlate to the intensity for the corresponding pixel in the other image. 
Other similarity measures including RMS and NCC have the expectation that the imaging modality used to generate both 
images are the same, and thus will not work when registering CT to MR images.
