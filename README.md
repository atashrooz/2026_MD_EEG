# Electrophysiological Correlates of Maladaptive Daydreaming

This repository contains the analysis code and derived datasets supporting the manuscript:

**Atashrooz, M., Kianimoghadam, A. S., Khosrowabadi, R., & Doosalivand, H.**  
*Electrophysiological Correlates of Maladaptive Daydreaming: An Exploratory Resting-State EEG Study*

The study examined resting-state EEG characteristics associated with probable maladaptive daydreaming (MD), with particular emphasis on spectral power, relative power, sensor-level functional connectivity, and frontal alpha asymmetry.

## Study overview

Resting-state EEG was recorded from 42 participants:

- 20 participants meeting the screening cutoff for probable maladaptive daydreaming
- 22 screen-negative participants

EEG was recorded during two resting-state conditions:

- Eyes open (EO)
- Eyes closed (EC)

Probable MD was defined using the Persian 16-item Maladaptive Daydreaming Scale (MDS-16), with a cutoff score of ≥ 50.

The main analyses focused on:

- Posterior alpha power
- Frontal/frontocentral theta power
- Frontal/frontocentral beta power
- Frontal theta/beta ratio
- Relative spectral power
- Sensor-level functional connectivity
- Frontal alpha asymmetry
- Exploratory electrode-level topographic analyses

## Repository contents

This repository includes MATLAB and R scripts used for EEG feature extraction, statistical analysis, and visualization, together with derived datasets used in the reported analyses.

Example structure:

```text
2026_EEG_MD/
│
├── README.md
│
├── data/
│   ├── eeg_band_power_features.csv
│   ├── eeg_connectivity_features_by_band.csv
│   └── eeg_participant_metadata.csv
│
├── preprocessing/
│   ├── create_connectivity_metadata.m
│   ├── datspec_extraction.m
│   ├── extract_connectivity_features_fieldtrip.m
│   ├── recompute_connectivity_features_6ROI.m
│   └── group_average_connectivity_visualization_v2.m
│
└── statistical_analysis/
    ├── spectral_regression_models.R
    ├── alpha_asymmetry_complete_analysis.R
    ├── eeg_connectivity_analysis.R
    └── visualize_connectivity_results.R
```

File names may differ slightly depending on the final repository organization.

## EEG acquisition and preprocessing

EEG was acquired using a 64-channel g.HIamp EEG system at a sampling rate of 512 Hz.

Preprocessing was performed in MATLAB using EEGLAB and included:

- Visual inspection of continuous EEG
- Band-pass filtering
- Identification and removal of persistently noisy channels
- Artifact Subspace Reconstruction (ASR)
- Independent component analysis (ICA)
- Removal of non-neural independent components
- Interpolation of excluded channels
- Common-average rereferencing

The scripts in this repository primarily reproduce the feature-extraction and statistical-analysis stages performed after EEG preprocessing.

## Spectral analysis

Power spectral density was estimated separately for eyes-open and eyes-closed recordings.

The following frequency bands were analyzed:

| Frequency band | Range |
|---|---|
| Delta | 1–<4 Hz |
| Theta | 4–<8 Hz |
| Alpha | 8–<13 Hz |
| Beta | 13–<30 Hz |
| Gamma | 30–45 Hz |

The principal regions of interest were:

### Posterior ROI

O1, Oz, O2, PO3, POz, PO4, P3, Pz, P4

### Frontal/frontocentral ROI

F7, F3, Fz, F4, F8, FC1, FCz, FC2, FC3, FC4

Absolute spectral-density values were converted to decibels prior to ROI-level statistical analysis.

Relative power was calculated from linear spectral-density values and analyzed as a complementary measure of the proportional distribution of EEG activity across frequency bands.

## Functional connectivity

Sensor-level functional connectivity was examined using:

- Magnitude-squared coherence
- Debiased weighted phase-lag index (wPLI)

Connectivity was estimated separately by resting-state condition and frequency band.

Regional connectivity analyses used six predefined scalp regions:

- Frontal
- Frontocentral
- Central
- Centroparietal
- Parietal
- Occipital

Within-region and between-region connectivity measures were calculated, yielding 21 regional connectivity measures in addition to a whole-scalp global measure.

Connectivity analyses were considered exploratory.

## Statistical analysis

Statistical analyses were conducted in R.

Primary spectral outcomes were analyzed using linear mixed-effects models including:

- Group
- Resting-state condition
- Group × condition interaction
- Participant-specific random intercept

The general form of the primary model was:

```text
EEG outcome ~ group * condition + (1 | participant)
```

False-discovery-rate correction using the Benjamini-Hochberg procedure was applied across predefined families of statistical tests.

Sensitivity analyses additionally examined whether the principal findings remained after adjustment for:

- Age
- Gender
- Visual imagery vividness
- Spontaneous imagery use
- Depression symptoms
- Anxiety symptoms
- Stress symptoms

## Software requirements

### MATLAB

The EEG analysis scripts require:

- MATLAB
- EEGLAB
- FieldTrip

EEGLAB version used in the study:

```text
EEGLAB 2025.1.0
```

### R

The R analysis scripts use packages including:

```r
tidyverse
readxl
lme4
lmerTest
emmeans
effectsize
broom
broom.mixed
psych
writexl
ggridges
```

Additional dependencies, where applicable, are indicated within individual scripts.

## Running the analyses

The general workflow is:

1. Obtain or prepare the preprocessed EEG recordings.
2. Run the MATLAB scripts for spectral feature extraction.
3. Run the MATLAB scripts for connectivity feature extraction.
4. Generate participant-level derived EEG datasets.
5. Run the R scripts for spectral statistical analyses.
6. Run the R scripts for alpha-asymmetry analyses.
7. Run the R scripts for functional-connectivity analyses.
8. Generate figures and supplementary visualizations.

The derived CSV files included in this repository allow the statistical analyses to be reproduced without rerunning the complete EEG preprocessing pipeline.

Local paths to EEG data, EEGLAB, and FieldTrip may need to be changed before running the MATLAB scripts.

## Data availability

The de-identified raw EEG data associated with this study are available on the Open Science Framework (OSF):

**OSF repository:** [ https://osf.io/gueqd ]
This GitHub repository contains the analysis code and derived datasets used for the spectral, connectivity, and statistical analyses reported in the manuscript.

The larger raw EEG files are hosted separately on OSF.

## Reproducibility

This repository is intended to provide a transparent record of the computational procedures used in the study and to facilitate reproduction of the reported analyses.

The derived datasets permit reproduction of the statistical models without requiring access to the original identifiable participant database.

Minor numerical differences may occur across versions of MATLAB, R, EEGLAB, FieldTrip, or individual R packages.

## Ethics

The study was conducted in accordance with the Declaration of Helsinki and was approved by the Research Ethics Committee of the School of Medicine, Shahid Beheshti University of Medical Sciences:

**IR.SBMU.MSP.REC.1402.512**

All participants provided informed consent prior to participation.

Only de-identified research data are intended for public sharing.

## Citation

If you use the code or data from this repository, please cite the associated manuscript:

> Atashrooz, M., Kianimoghadam, A. S., Khosrowabadi, R., & Doosalivand, H.  
> *Electrophysiological Correlates of Maladaptive Daydreaming: An Exploratory Resting-State EEG Study.*

Publication details and DOI will be added following publication.



## Contact

For questions regarding the study or repository:

**Amir Sam Kianimoghadam**  
Department of Clinical Psychology  
School of Medicine  
Shahid Beheshti University of Medical Sciences  
Tehran, Iran

Email: as.kianimoghadam@gmail.com
