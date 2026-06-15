# Check Excess Params

Function for checking excessive parameter names. This function will
notify the user if any unavailable parameters are called for any given
given function. This function is designed to work as internal
function-call in already available GAMBLR functions.

## Usage

``` r
check_excess_params(...)
```

## Arguments

- ...:

  Parameters to check.

## Value

Nothing

## Details

Catch function calls containing unsupported arguments.
