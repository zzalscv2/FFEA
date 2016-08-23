FFEA Analysis {#FFEAanalysistut}
=============================

## The FFEA python package

In previous parts of this tutorial, we have employed several python scripts by invoking 'ffeatools' at the terminal. These scripts all make use of a set of core FFEA python modules. Each module corresponds to a different type of data - the trajectory, pin files, the stokes file, et cetera. These are all unified and linked by the FFEA script file.

In addition to a few basic analysis tools, these core Python modules are provided as part of a Python package, which allows them to be imported into a Python interpreter. To install this package, open the FFEA source folder (containing setup.py) in the terminal, and type

```sh
python setup.py install
```

This will install the FFEA tools into your python site-packages folder. For this tutorial, you may want to consider installing the [Anaconda Python distribution](https://www.continuum.io/downloads), which comes with a set of common scientific python packages, and allows for packages to be installed and removed easily. An interactive python shell (such as IPython) or a Python IDE with an object inspector and auto-complete is recommended. Anaconda comes with one called Spyder, which can be launched using the terminal

```sh
spyder
```

In our interactive session we first
```python
import ffeatools
```
The FFEA tools can now be accessed in the ffeatools namespace. To start, we create a new script object.

```python
our_script = ffeatools.modules.script("/path/to/script.ffea")
```

We can now load in the contents of any of the files that the script is linked to. We will start by loading the trajectory.

```python
our_trajectory = our_script.load_trajectory()
Loading FFEA trajectory file...

Frames read = 1, Frames skipped = 0
Frames read = 2, Frames skipped = 0
Frames read = 3, Frames skipped = 0
...
Frames read = 99, Frames skipped = 0
Frames read = 100, Frames skipped = 0
Frames read = 101, Frames skipped = 0
done! Successfully read 101 frame/s from '/path/to/trajectory.out'.
```
We can now access all of the data inside the trajectory. For example,

```python
our_trajectory.num_blobs
1
```python
