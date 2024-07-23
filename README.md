# SWIFTJ1727.8-1613 QPO Analysis

I am analyzing light curves from the observations of SWIFTJ1727.8-1613 for Quasi Periodic Oscillations in an indirect analysis of the system and its properties.

## Installations and Getting Started

Install Anaconda and create a new environment for this project

```
conda create -n <env name> batanalysis matplotlib astropy numpy scipy boto3 sgp4
```

Then to enter:

```
conda activate <env name>
```

### BatAnalysis Installation

```
pip install batanalysis
```

```
git clone https://github.com/parsotat/BatAnalysis.git
cd BatAnalysis/
pip install -e .
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

To run notebook, make sure Jupyter Notebook is installed and run
```
jupyter notebook J1727.ipynb
```


## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE.md file for details