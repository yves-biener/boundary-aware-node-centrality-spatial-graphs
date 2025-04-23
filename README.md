# Boundary-aware node centrality spatial graphs

University project for boundary-aware node centralities for spatial graphs. For the project details and results see the [typst paper](./doc/paper.typ). For building the document and the usage of the python library used for this project see the corresponding sections below.

## Code usage

> [!IMPORTANT]
> The implementation uses [gurobi](https://www.gurobi.com/solutions/gurobi-optimizer/) for solving the linear problem for the function fitting which requires a *license*. Please refer to the gurobi license documentation for details.

Install the requirements into a virtual environment:

```sh
# create virtual environment
python -m venv venv
# activate virtual environment
source venv/bin/activate
# install required dependencies
pip install -r requirements.txt
```

The `src` directory contains all the python source files that are used for the creation of the diagrams and images found in the `doc/figures` directory. For the corresponding usages for both datasets used for this project please see [merfish.py](./merfish.py) and [mibitof.py](./mibitof.py).

For running the generation of the diagrams and images you can run one of the two python scripts like so (inside of the virtual environment):

```sh
python merfish.py
```


## Document generation

The document is written using [typst](https://typst.app) and can be compiled into a pdf file using:

```sh
typst c doc/paper.typ
```

The generated pdf file will then be generated into the `doc` directory named `paper.pdf`.
