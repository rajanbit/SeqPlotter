# SeqPlotter: Python package for sequence data analysis and visualization

### Introduction
SeqPlotter is a Python package developed to streamline the visualization and analysis of biological sequence data. It provides
seamlessly integration into your Python workflow, with robust data handling capabilities, intuitive visualization, and useful
analysis features. From parsing raw sequence files to creating stunning plots and charts.

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
pip install dist/SeqPlotter-0.1.0-py3-none-any.whl
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
