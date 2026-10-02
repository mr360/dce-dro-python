# dce-dro-python
 Python version of a DCE DRO phantom. Based on the matlab version available at https://github.com/anstepsa/dce-dro-matlab

# Virtual environment
I use ```mkvirtualenv``` on MacOSx (Sonoma 14.5). To set it up, instructions can be found elsewhere, such as in the [official documentation](https://virtualenvwrapper.readthedocs.io/en/latest/index.html). I haven't tested other options, but I believe there is no much difference when using other tools, such as the python native [```venv```](https://docs.python.org/3/library/venv.html). 

Once you are in the virtual environment, install the required modules listed in the ```requirements.txt```:

```
pip install -r requirements.txt
```

# FTV validation
The [`validation/`](validation/README.md) folder contains tools to validate and audit the FTV, PE and SER values computed by the SeQ-DCEMRI 3D Slicer extension with the FTV digital phantom. See [`validation/METHODOLOGY.md`](validation/METHODOLOGY.md).
