# Change gains on the loss

Dec 9/2022
current loss =  doubleHinge_x  + 12.0f0*(1.0f0-cos(x2)) + 2.0f0x1dot^2.0f0 + 0.5f0*x2dot^2.0f0
TODO: - the weights on the loss are not balanaced.
With these gains, the cart does not want to swing but the pendulum stabilizes
when close to the upright.
- Try lowering the cost on x1dot to allow it to swing up.
Moreover, the control NN may be too big now. It has large variance in control across the 
states.

Jan 12/2023
current loss = 5.0f0*(doubleHinge_x  + 8.0f0*(1.0f0-cos(x2)) + 1.0f0*abs(x1dot) + 0.5f0*abs(x2dot))

This allows it to pump. Catching against the wall is not shown yet. The positive x2dot space
seems to be filled with only one bin. This does not seem rights

Feb 6/2023
Set distance loss has much easier time pumping. Catching is not as robust.
Set distance loss is min(8.0f0*(1.0f0-cos(x2)) + 1.0f0*abs(x1dot) + 0.8f0*abs(x2dot)) + doubleHingeLoss.

-Ever since adding the damping and the low restitution, it has not needed a third bin.

