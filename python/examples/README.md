There are currently three examples in this directory.

(1) `run_simple_example.py` is the complete code for a very simple example which applies Anderson Acceleration to a small problem (defined within the file). To run:

```python
python run_simple_example.py
```

(2)`run_poisson_example.py` is the driver script tp apply AA to a nonlinear Poisson equation. The fixed point function is defined in `timestepperfunc_poisson.py`. This example can be run in serial or parallel. For serial:

```python
python run_poisson_example.py
```

For parallel, modify the script so that `useMPI=True` (toward the top of the file) and do:

```python
mpiexec -n 2 python run_poisson_example.py
```

(3) The third example (run_poisson_xy_example.py) solves the same problem as the second, except that it splits the state vector x into two vectors, one confusingly also called x and a second called y. The AA extrapolation coefficients are calculated using x and the corresponding residual, but also applied to the y vector. The fixed point function is defined in `timestepperfunc_poisson_xy.py`. You can run this example in serial or parallel similarly to the previous example.
