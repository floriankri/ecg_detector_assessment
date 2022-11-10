# Assessment of Electrocardiogram Beat Detectors Using Synthetic and Real-World Data

This repository contains three main resources. **Algorithm_tester** contains databases, algorithms and a script to test them. **Annotator** has a MATLAB script that allows annotation of WFDB ECG signals. **Signal_generator** contains MATLAB scripts that can generate synthetic ECG signals.

## Authors
- Florian Kristof [floriankri](https://github.com/floriankri)
- Peter Charlton [peterhcharlton](https://github.com/peterhcharlton)
- Leon Nissen
- Maximilian Kapsecker

## Citation

## Installation

### Algorithm Tester
1. Install **Visual Studio Code**
2. Install **MATLAB R2022a** or higher
3. Install **MATLAB Engine API for Python** available [here](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html) in Visual Studio Code
4. `pip install requirements.txt` in Visual Studio Code
5. Try executing `/Algorithm_tester/main.ipynb` and see if and where it runs into errors

### Annotator
1. Install **MATLAB R2022a** or higher
2. Install **Waveform Database Software Package (WFDB) for MATLAB and Octave** available [here](https://physionet.org/content/wfdb-matlab/0.10.0/)
3. Change all permanent paths in the files to your local paths
4. Change paths to signals that should be used in `/Annotator/MAIN.m`

### Signal Generator
This part is automatically executed by **Algorithm Tester**. If you want to setup the unchanged version please visit the initial page called [Model for Simulating ECG and PPG Signals with Arrhythmia Episodes](https://physionet.org/content/ecg-ppg-simulator-arrhythmia/1.3.1/).

## Usage and Results

### Algorithm Tester

#### Comparing Atrial Fibrillation and Sinus Rhythm
![Comparing Atrial Fibrillation and Sinus Rhythm](./Algorithm_tester/figures/af_sr_comparison_v001.svg)

### Annotator

### Signal Generator