# SchrodingerEq_1D_tutorial
Solving Schrodinger Equation Numerically in Python

## Files
`SchrodingerEq_1D_tutorial.ipynb`  A tutorial of solving 1-D Schrodinger equation. (Simple-Harmonic-Oscillator)
`step_potential/step_potential.py` A script to generate the animation that shows a wave-packet scattered by a step-potential.

## Environment (macOS, suggested)
1. Install [MacPort](https://www.macports.org/install.php)
2. Install **ImageMagick** (for GIF) and **ffmpeg** (for MP4)
```
sudo port install ImageMagick ffmpeg +gpl +postproc +lame +theora +libogg +vorbis +xvid +x264 +a52 +faac +faad +dts +nonfree
```
3. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
4. Install Python packages
```
conda install numpy scipy matplotlib jupyter ipython
```
5. `conda activate`
6. `jupyter-notebook`
