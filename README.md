# 10-minimizers

Accompanying code for the paper *10-minimizers: a promising class of constant-space minimizers*.

## Project Structure
### `density-measurement-and-plots/`
This directory contains the interactive Python notebooks used to benchmark the density of various minimizer schemes.

* **`minimizer_density_comparison.ipynb`**: The primary notebook for running benchmarks and generating plots.
    * The first cell specifies the library versions and Python environment used.
    * **Compatibility**: Designed to run on most modern Python versions with Jupyter (`.ipynb`) support.
* **`greedymini_density.ipynb`**: Specifically used to measure the density of **GreedyMini** minimizers.
* **`runtime_plots.ipynb`**: Used to generate key-retrieval graphs. 
    * *Note: The core calculations are performed in C++; this notebook handles the visualization of those results.*
### `key-retrieval-code/`
Contains the C++ implementation used to measure key retrieval performance across different minimizer schemes.
> [!NOTE]
> While functional for benchmarking, there is room for optimization in the implementation, particularly for the baseline minimizers we did not develop.

## Requirements

### Python
To open and run the interactive Python notebooks (`.ipynb`), you will need one of the following:

* **Jupyter Notebook or JupyterLab**: The standard web-based environments for interactive computing.
* **VS Code**: With the **Jupyter** extension installed.
* **Google Colab**: A cloud-based option that requires no local setup.

#### Python Environment
* **Python Version**: Developed and tested using **Python 3.11.4**.
* **Core Libraries**:
    * `numpy`: 2.4.2 
    * `pandas`: 2.3.0 
    * `matplotlib`: 3.10.3 
    * `tqdm`: 4.67.1 
    * `joblib`: 1.3.2
> [!NOTE]
> While these specific versions were used for benchmarking, we expect the code to remain functional across most modern versions of these libraries that are compatible with one another.
#### Quick Setup
You can install all necessary dependencies using `pip`:

```bash
pip install numpy==2.4.2 pandas==2.3.0 matplotlib==3.10.3 tqdm==4.67.1 joblib==1.3.2 notebook
```

### C++
To compile the code in `key-retrieval-code/`, ensure you have:
* **C++ Compiler:** Must support **C++20** or later (.e.g, `g++-10+` or `clang-10+`)
* **Build Command:** When compiling manually, ensure you use the `-std=c++20` flag.

#### Example Compilation
```bash
g++ -O3 -std=c++20 Minimizers.cpp sampling_runtimes.cpp -o key_retrieval
```

## Contact
In case of issues, you may contact us at tziony.i@gmail.com.
