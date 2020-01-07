# Interactive 2D Gravity Simulation in Javascript
https://jurasofish.github.io/gravity/

Click and drag to place a new body.

Use spacebar to pause and right arrow to step through time.

Timescale and granularity of the simulation can be controlled with the &Delta;t and Lookahead inputs.

Try setting the Follow dropdown to "Earth" to see the Moon around the Earth.

Scroll to zoom.

This is implemented by solving the n-body problem (Newton's law of gravitation) with [Ricky Reusser's javascript ODE45 solver](https://github.com/scijs/ode45-cash-karp).

Collisions are perfectly inelastic.

Careful it can hang the browser if you make the system do something tricky :)
