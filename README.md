A tool to compute the exact minimum circuit-level distance of a QEC code given in stim circuit format using integer linear programming.

A detailed explanation is given in this blog [post](https://shohamjac.github.io/posts/circuit-Level-distance/).

example:
```python
python BB_code_example.py
```

prints:
```
Solver: CBC
 | distance:  12
 | optimal:  True
 | time:  78.36955916206352

Solver: GRB
 | distance:  12
 | optimal:  True
 | time:  3.3316898189950734
```