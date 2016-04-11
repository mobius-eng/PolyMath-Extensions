I am Conjugate Gradient method with EWT-RMS residual error: the residual error is calculated as:

``` 
      __________________________
     / /          r           \  2
    / 1|           i          | 
|  /  -|----------------------| 
| /   N|RelTol  v   +  AbsTol | 
|/     \      i  i           i/ 


```

where `r_i` are the components of the residual, `RelTol` and `AbsTol` are relative and absolute tolerance vectors and `v` is a reference vector (by default equal to RHS vector).

