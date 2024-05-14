[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10645748.svg)](https://doi.org/10.5281/zenodo.10645748) 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



# SeqPlotter: Python package for sequence data analysis and visualization

### Introduction
SeqPlotter is a Python package developed to streamline the visualization and analysis of biological sequence data. It provides
seamlessly integration into your Python workflow, with robust data handling capabilities, intuitive visualization, and useful
analysis features. From parsing raw sequence files to creating stunning plots and charts.

<p align="center" width="100%"> <img width="90%" src="https://raw.githubusercontent.com/rajanbit/SeqPlotter/main/docs/img0.png"></p>

### Requirements
- Python >=3.8
- Git

### Installation
**pip**
```
git clone https://github.com/rajanbit/SeqPlotter.git
cd SeqPlotter/
python -m pip install --upgrade build
python -m build
pip install dist/SeqPlotter-0.3.0-py3-none-any.whl
```

**conda**
```
conda create -n seqplotter python=3.8 git
conda activate seqplotter
pip install git+https://github.com/rajanbit/SeqPlotter.git#egg=SeqPlotter
```

### Documentation
- User Guide: https://github.com/rajanbit/SeqPlotter/blob/main/docs/SeqPlotter_user_guide.pdf
- Tutorials: https://github.com/rajanbit/SeqPlotter/blob/main/docs/SeqPlotter_tutorial.ipynb

### Citation
Rajan. (2024). SeqPlotter: Python package for sequence data analysis and visualization (v0.1.0). Zenodo. https://doi.org/10.5281/zenodo.10645748
