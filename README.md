# BBAND Index and Adaptive-Debanding Filter

This repository contains reference matlab implementations of the BBAND Index
and Adaptive-Debanding Filter from the following papers:

1. Z. Tu, J. Lin, Y. Wang, B. Adsumilli and A. C. Bovik, "BBAND INDEX: A NO-REFERENCE BANDING ARTIFACT PREDICTOR," ICASSP 2020 - 2020 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), Barcelona, Spain, 2020, pp. 2712-2716, doi: 10.1109/ICASSP40776.2020.9053634 [IEEE Explore](https://ieeexplore.ieee.org/document/9053634).
1. Z. Tu, J. Lin, Y. Wang, B. Adsumilli, and A. C. Bovik. "Adaptive Debanding Filter". IEEE Signal Processing Letters, Volume 27. [IEEE Explore](https://ieeexplore.ieee.org/document/9201351).

# Running the code

There are third_party routines in "third_party/matlabfns". You will need to add
that directory to your path (e.g., path(path(), "${PATH_TO_ROOT}/third_party/matlabfns")).

The code for the blind banding detector (BBAND) and the adaptive debanding filter
are in the following files.
```
src/BBAD_I.m
src/deband_filter_yuv420p
```

Please see the comments in the files for description of the input/output parameters.

Alternatively, you can have a look a scripts used to run experiments in the papers:

1. The script for running banding detection experiments is located in src/video_band_detect.m.
1. Debanding filtering as preprocessing experiments are in src/video_preproc_dither.m and
the script for running adaptive debanding postprocessing is src/video_postproc_deband_filter.m.



## Tests

There are also a few tests in the "tests" folder that illustrate how to call the routines.
You can run the tests from the  by adding the "src" directory to your path.

```
test_BBAD_I
test_deband_filter
```

The scripts will output whether the implementation passed or failed. E.g., at the end
of the debanding script you should see lines like:

```
KristenAndSara_10frames_512x256 passed
```







