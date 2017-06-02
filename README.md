# Simple IPOPT wrapper using finite differences

This project uses ipopt for a simple optimization problem using finite differences to calculate the objective function gradient and the constraint jacobian.

[The optimization test case](http://ab-initio.mit.edu/wiki/index.php/NLopt_Tutorial "NLP Example")

## Tested on:

- Windows 10 - mingw64

## How to use:

- Create a class inhering from `OptInterface` and define some optimization parameters and functions as x low and upper bounds; gradient low and upper bounds; objective function (override) and constraint function (override). In this case it is the `optexample1.h/myExample1`.

- Create the `main` function, which defines the case example, the Ipopt C++ wrapper class `TNLP`. To run the optimization call `createAppAndRun(nlp);`

## TODO:

- Test on debian
- Check how the Ipopt was installed

*(Note: was used the [newcppproject.sh](https://gist.github.com/caiofcm/83d4d3d2370546d846454ff74dea7348) to generate the Makefile)*
