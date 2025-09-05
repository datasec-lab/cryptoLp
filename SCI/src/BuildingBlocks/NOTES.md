#### Findings on AuxProtocols::multiplexer

**Objective.**


```plaintext
y = [sel_0*x_0 + sel_1*(x_0 - 2*sel_0*x_0)] + [sel_1*x_1 + sel_0*(x_1 - 2*sel_1*x_1)]
```

**Step 1.** Each party computes `corr_data_b = x_b - 2*sel_b*x_b`.

**Step 2.** Party `b` sends `corr_data_b` to party `1-b`, and it holds `data_S_b`.

**Step 3.** Party `b` receives `data_R_{1-b}` from party `1-b`.

**Step 4.** Each party computes `y_b = sel_b*x_b - data_R_{1-b} - data_S_b`.


Note that `y_0 + y_1` equals to the objective output.


**TODO.**
- Figure out the purpose of `iknp_reversed` and `iknp_straight`.