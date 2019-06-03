# ECG signal processing

![ecg](https://github.com/trevorwitter/ECG-processing/blob/master/images/ecg.png)

Tools for processing raw electrocardiogram data

#### QRS Detection
The QRS complex is the most prominent feature in ECG signal and is therefore useful for quantifying timing of individual heartbeats. This allows for quantifying useful physiological metrics such as R-R interval, heart rate, heart rate variability as well as providing timing for measurment of other signals such as beat-to-beat systolic/diastolic blood pressure and brain blood flow velocity.


The `ecg_wavelet()` class in ecg_processing.py provides a [butter lowpass filter](https://en.wikipedia.org/wiki/Butterworth_filter) for noise removal and QRS detection via [wavelet transform](https://en.wikipedia.org/wiki/Wavelet_transform). The `.get_qrs()` method returns a dataframe including a column labeling each detected QRS timepoint.  

![QRS detect](https://github.com/trevorwitter/ECG-processing/blob/master/images/qrs_detect.png)
