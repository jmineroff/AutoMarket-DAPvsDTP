# AutoMarket-DAPvsDTP
Implementation of Design-and-Pricing (DAP) and Design-then-Pricing (DTP) portfolio strategies for the automotive market using Matlab.

Corporate profits in a simple 2 company market (with one vehicle each) are modeled under both portfolio strategies. In the DAP model, manufacturers design a product portfolio without planning for the effect of future price competition on profits; DTP adds this consideration. The performance of a direct formulation and a more robust formulation are compared with multiple optimization solvers.

## Dependencies
This program uses solvers from the [Matlab Optimization Toolbox](https://www.mathworks.com/products/optimization.html).

## File Organization
* 'OPT....m': Contain the high level optimization code to compare the effects of the Zeta map on interior point or sequential quadratic programming solvers
* 'OBJ...m': Objective function definitions for different problem formulations
* 'NONLINCON...m': Nonlinear constraint function definitions

## Publications
1. W Ross Morrow, Joshua Mineroff, Kate Whitefoot (2014). Numerically Stable Design Optimization with Price Competition. _Journal of Mechanical Design_
1. W Ross Morrow, Joshua Mineroff, Kate Whitefoot (2012). Numerically Stable Design Optimization with Price Competition. _American Society of Mechanical Engineers: International Design and Technical Conference_
3. Kate Whitefoot, Meredith Fowlie, Steven Skerlos (2011). Product Design Responses to Industrial Policy: Evaluating Fuel Economy Standards using an Engineering Model of Endogenous Product Design. _Energy Institute at Haas Working Paper Series_

## Licensing
This code is licensed under MIT license.
