CISC 472: Assignment 3
Evelyn Yach (20071956)
Daniel Oh (20063998)

Results(Without RANSAC):
pivot_calibration_0: (-3.62744, -10.8191, -186.279) +/- (16.6282, 14.1962, 12.5941)
pivot_calibration_1: (-2.95554, -18.0382, -104.901) +/- (41.7893, 31.2551, 20.855)
pivot_calibration_2: (-17.3991, -0.786382, -158.162) +/- (0.449139, 0.309227, 0.143006)
pivot_calibration_3: (-17.9292, -0.639105, -160.257) +/- (0.811305, 0.614405, 0.107862)

Results(With RANSAC):
pivot_calibration_0: (-18.1271, -1.04205, -159.227) +/- (1.44107, 1.73756, 0.327935)
pivot_calibration_1: (-17.1314, -0.902285, -157.934) +/- (0.43045, 0.427003, 0.190469)
pivot_calibration_2: (-17.3961, -0.970535, -158.61) +/- (0.669941, 0.539175, 0.183115)
pivot_calibration_3: (-17.9205, -1.01297, -160.639) +/- (0.85923, 0.690126, 0.139692)


5a. What characterizes a good calibration set?  Justify based on your results.

A good calibration set would be one where there are minimal outliers from the surface we’re 
calibrating to. This is demonstrated as RANSAC which removes those outliers and performs better 
with the data. If our number of outliers is minimized then we don’t have to worry about throwing 
our calculations off to begin with.

   
5b. What method (sphere or sphere-with-ransac) was best for handling each calibration set?  
Is it always the same method, or do different characteristics of the set deserve different methods? 
Justify your answer.

In pivot_calibration_0 and pivot_calibration_1, the sphere-with-ransac method is clearly better. 
However, for pivot_calibration_2 and pivot_calibration_3, the difference becomes more ambiguous. 
Looking at standard deviation, we can see that the non RANSAC method actually performs better in 
pivot_calibration_2. This could be due to the fact that there are almost no egregious outliers in 
pivot_calibration_2. In the RANSAC method, our code is removing some points which may actually be 
okay to include. Using those additional points, it’s possible that we are achieving a tighter 
standard deviation.

In conclusion, the sphere-with-ransac method can be used in every scenario and get a sufficiently 
accurate result. The sphere (no RANSAC) method should only be used in scenarios where there are no 
outliers in the dataset.


5c. What recommendations do you have for pivot calibration?

Starting the pivot calibration process when the tip of the tool is already in the swivel point 
where you plan on calibrating the tool is a good way to avoid outliers. Another good practice would 
be to try and achieve high coverage on possible inliers by starting with swiveling the tool in 
large circles, and gradually shrinking the size of the circles as the process continues. This would 
create a spiral pattern of pivot calibration points that can cover a large calibration area.