# SWIFTJ1727.8-1613 QPO Analysis

I am analyzing light curves from the observations of SWIFTJ1727.8-1613 for Quasi Periodic Oscillations in an indirect analysis of the system and its properties.

## Installations and Getting Started

Make sure you have installed  Anaconda and create a new environment for this project

```
conda create -n astro batanalysis matplotlib astropy numpy scipy boto3 sgp4 pathlib stingray lmfit
```

Then to enter:

```
conda activate astro
```

## Running

```
python3 swiftJ1727_lightcurve.py
```

For simple visualization

```
python3 plot_counts.py
```

A Fourier Transform of the light curve
```
python3 Fourier.py
```

To run the notebook, first install Jupyter Notebook by running 
```
python -m pip install jupyter
```

Then to open the Jupyter notebook for this project run:

```
jupyter notebook J1727.ipynb
```

## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details