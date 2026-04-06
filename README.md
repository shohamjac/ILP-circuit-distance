A tool to compute the exact minimum circuit-level distance of a QEC code given in stim circuit format using integer linear programming.

A detailed explanation is given in this blog [post](https://shohamjac.github.io/posts/circuit-Level-distance/).

example:
```bash
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

### Calculating the distance of BB18 from 2603.28627

In 2603.28627, they use a BB code with the parameters $[[248,10,\le18]]$.

Using those scripts, we can easily check that 18 is indeed tight. Use:
```bash
python oratomic_BB.py --basis z --time-limit 10000
```

It should run for about half an hour (on a 32-cores machine) and return:
```
Solver: GRB
 | distance:  18
 | optimal:  True
 | time:  1572.8003788780188
```

Do it for both bases to get both $d_x$ and $d_z$.