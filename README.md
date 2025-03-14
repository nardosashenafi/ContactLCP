# ContactLCP
This repository contains a Julia package that resolves contacts via Linear Complementarity Problem (LCP). Moreover, the package allows auto-differentiation through the LCP solver. Examples are provided on how to use auto-differentiation capabilities to estimate system parameters and learn optimal controllers through data-driven techniques.

The package includes:
- LCP solver (Lemke's algorithm)
- Models of the rimless-wheel, bouncing ball and cart-pole with soft walls
- System parameter estimation for the bouncing ball
- Data-driven scheme to learn control of the rimless-wheel to achieve constant hip speed
- Bayesian learning technique to learn switching controller for the cart-pole with softwalls

Rimless wheel: controlled via data-driven passivity-based control (NeuralPBC)
![ezgif com-gif-maker](https://user-images.githubusercontent.com/36546726/197265071-11888428-06d2-49e2-b9de-3aca88f47fed.gif)

Cart-pole with soft wall contact: controlled via LQR

![ezgif com-gif-maker(1)](https://user-images.githubusercontent.com/36546726/203182031-972f14d3-be52-440b-9b0e-7909e83a93ab.gif)

Cart-pole swing up with Deep-net Mixture of Experts

![MOESimExpPlots](https://user-images.githubusercontent.com/36546726/218184094-8ead79bc-fce8-4901-b2ca-6fc32d9bdebb.png)
