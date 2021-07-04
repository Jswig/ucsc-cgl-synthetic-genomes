# Mixture model

## Symmetric Mixture Model 

### Log approximation

$$
\log(x + y) = \log(\exp(\log x) + \exp(\log y))
$$
$$ 
          \approx \log(\exp(\log x) * (1 + \exp(\log y - \log x)))
$$
$$
           \approx \log x + \log(1 + \exp(\log y - \log x))
$$
and then you generally make sure x is the larger of the two variables (you end up sacrificing precision on y, so you want to make sure it’s the less significant variable)

Can apply repeatedly with x as the cumulative sum

*Sum of arbitary number of terms*

$\log(\sum_{i=1}^N x_i) = \log(\sum_{i=1}^{N-1} x_i) + \log(1 + \exp(log x_n - \log \sum_{i=1}^{N-1} x_i))$

In python 

```python
def logsum_approx(summands):
    n = len(summands)
    sum = 0
    while (n > 0):
            
        n -= 1
```


https://github.com/benedictpaten/sonLib/blob/master/C/impl/sonLibMath.c