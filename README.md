# SchrodingerEq_1D_tutorial
Solving Schrodinger Equation Numerically in Python.  
Please refer to [this blog article](https://blog.gwlab.page/solving-1-d-schrodinger-equation-in-python-dcb3518ce454) for more detail explaination of the code.

## Files
- `SchrodingerEq_1D_tutorial.ipynb`  A tutorial of solving 1-D Schrodinger equation. (Simple-Harmonic-Oscillator)
- `step_potential/step_potential.py` A script to generate the animation that shows a wave-packet scattered by a step-potential.
- `step_potential_cpp/step_potential_rk4.cpp` The C++ source code for simulating a wave-packet scattered by a step-potential.
- `step_potential_cpp/plot.py` A script to generate the animation from the result of step_potential_rk4.

## Environment (macOS, suggested)
1. Install [MacPort](https://www.macports.org/install.php)
2. Install **ffmpeg** (for MP4)
```
sudo port install ffmpeg +gpl +postproc +lame +theora +libogg +vorbis +xvid +x264 +a52 +faac +faad +dts +nonfree
```
3. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
4. Install Python packages
```
conda install numpy scipy matplotlib jupyter ipython imagemagick
```
5. `conda activate`
6. `jupyter-notebook`

## Environment for C++ version (macOS, suggested)
1. Install [Eigen](https://eigen.tuxfamily.org/)
```
# macOS
sudo port install eigen3

# Ubuntu
sudo apt install libeigen3-dev
```
2. Modify `step_potential_cpp/Makefile` with proper path to your Eigen installation.
3. Compile
```
cd step_potential_cpp
make
```
3. Run and generate plots(animations)
```
cd step_potential_cpp
./step_potential_rk4
conda activate
./plot.py
```
