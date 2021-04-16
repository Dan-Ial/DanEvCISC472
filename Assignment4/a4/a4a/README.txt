CISC 472: Assignment 4a
Evelyn Yach (20071956)
Daniel Oh (20063998)

2.
knee1.csv

    Rotation:
        [ 0.6991   -0.6417    0.3155;
          0.1763    0.5822    0.7937;
         -0.6930   -0.4992    0.5202 ]
    Translation:
        [ 616.6926  823.0078  546.7715 ]
    Num Of Iterations:
        46
    RMS Error:
        2.8958

knee2.csv

    Rotation:
        [ 0.8336    0.1136    0.5405;
          0.3191    0.6996   -0.6393;
         -0.4508    0.7054    0.5470 ]
    Translation:
        1.0e+03 * [ 0.7675   -1.3357    0.2858 ]
    Num Of Iterations:
        42
    RMS Error:
        2.0853


4.
knee1.csv

    Num Of Attempts:
        10
    Rotation:
        [ -0.2937    0.9328   -0.2090;
          -0.9195   -0.2159    0.3284;
           0.2612    0.2887    0.9211 ]
    Translation:
        1.0e+03 * [ -505.4587  333.2505  948.8594 ]
    Num Of Iterations(of best attempt):
        51
    RMS Error:
        1.5684

knee2.csv

    Num Of Attempts:
        10
    Rotation:
        [ 0.6996    0.7021    0.1323;
         -0.2897    0.4480   -0.8458;
         -0.6531    0.5534    0.5169 ]
    Translation:
        1.0e+03 * [ 0.0446   -1.5753    0.2778 ]
    Num Of Iterations(of best attempt):
        31
    RMS Error:
        1.8388


5. Did the multi-attempt version of ICP work better?  Discuss this
briefly, speculating on ICP's strengths and weaknesses (at most a
paragraph).

The multi-attempt version of ICP consistently produced better RMSE values
than the single attempt version. With single attempt, there is a chance that
the initial translation and rotation can lead to an undesirable local minimum.
With multi-attempt, starting with randomly generated initial translations and 
rotations can eliminate this downside. However, runtime is much longer 
depending on the number of attempts the user desires. With a low number of 
attempts, we can possibly run into the issue of achieving an undesirable local
minimum.






