# Learning Switching Control for the cart-pole with soft wall contact

* The objective of this work is to learn three controllers that stabilize the 
cart-pole bounded by two walls **D** units apart. The contact from the walls create
disturbances that cannot be trivially handled by classic controllers such as LQR. 
The switching controller allows the system to respond differently with the 
discontinuities in the system dynamics. 

**The Drunken cart-pole** : two neural network controllers and a gating network inferred via expectation maximization

![ezgif com-gif-maker(4)](https://user-images.githubusercontent.com/36546726/204888204-78ee491c-4098-4803-b389-630164a51e33.gif)

The heat map on the top left shows that when the angle of the pendulum is outside the range of [0, 2] radians, it applies the first controller that swings up the pendulum. When the pendulum angle is within [0, 2] radians and it hits the wall, the second controller activates and balances the pendulum agains the wall. 

**Cart-pole switching between three controllers** :

![ezgif com-gif-maker(5)](https://user-images.githubusercontent.com/36546726/205130570-09daa14e-28a0-45bc-9892-7eee0de9032a.gif)



