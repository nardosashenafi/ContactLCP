# Change gains on the loss

current loss =  doubleHinge_x  + 12.0f0*(1.0f0-cos(x2)) + 2.0f0x1dot^2.0f0 + 0.5f0*x2dot^2.0f0
TODO: - the weights on the loss are not balanaced.
With these gains, the cart does not want to swing but the pendulum stabilizes
when close to the upright.
- Try lowering the cost on x1dot to allow it to swing up.
Moreover, the control NN may be too big now. It has large variance in control across the 
states.

