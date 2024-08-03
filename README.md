# SWIFTJ1727.8-1613 QPO Analysis

I am analyzing light curves from the observations of SWIFTJ1727.8-1613 for Quasi Periodic Oscillations in an indirect analysis of the system and its properties.

## Installations and Getting Started
First clone this repository with:

```
git clone git@github.com:Dirac12/SWIFT-Quasi_Periodic_Oscillations.git
```


First install Anaconda, where you can find instructions here: <a href="https://docs.anaconda.com/anaconda/install/linux/" target="_blank">InstallAnaconda</a>

It is additionally recommended that you create a new environment for this project. Run the following in a terminal

```
conda create -n astro batanalysis matplotlib astropy numpy scipy boto3 sgp4 pathlib stingray lmfit
```

Then to enter:

```
conda activate astro
```

## Running

Enter the SWIFT directory by running
```
cd SWIFT-Quasi_Periodic_Oscillations
```

To download initial data run the following in a terminal:
```
python3 swiftJ1727_lightcurve.py
```

For simple visualization

```
python3 plot_counts.py
```

And to view a Fourier Transform of the light curve
```
python3 Fourier.py
```

To run the notebook, first install Jupyter Notebook by running 
```
python -m pip install jupyter
```

Then to open the Jupyter Notebook for further analysis
for this project run:

```
jupyter notebook J1727_4.ipynb
```

## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details