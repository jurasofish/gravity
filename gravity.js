"use strict";

const G = 6.67408e-11; // gravitational constant
const TOL = 1e-8;  // Tolerance for ODE solver.
let canvas = document.getElementById("myCanvas");
let ctx = canvas.getContext("2d");
let HEIGHT = 400;
let WIDTH = 600;
canvas.height = HEIGHT;
canvas.width = WIDTH;


let drawLine = function(ctx, pts) {
    ctx.beginPath();
    pts.forEach(point => {
        ctx.lineTo(point[0], point[1]);
    });
    ctx.stroke();
    ctx.closePath();
}

let col_2d = function(x, n) {
    /* given a 2D array x, return column n as a 1D array */
    return x.map(row => row[n]);
}

let max_1d = function(x) {
    /* given a 1D array x return the maximum value */
    return x.reduce((a, b) => Math.max(a, b));
}

let min_1d = function(x) {
    /* given a 1D array x return the minimum value */
    return x.reduce((a, b) => Math.min(a, b));
}

let Body = class{
    /* A body is a physical object which produces and is affected by gravity. */
    constructor(name, m, p, v, a, r, t=0) {
        this.name = name; // Name
        this.m = m; // mass (kg)
        this.a = a; // array of x, y of applied acceleration (m/s/s)
        this.r = r; // Body radius, for for collissions (m)
        this.r_g = Math.log(r) * 5e8; // Body radius, for graphics.
        
        // Store the expected future position over time.
        // Similar to t_hist and p_hist.
        this.t = [t];  // 1D array of time (s)
        this.p = [p];  // 2D array of [x, y] (m)
        this.v = [v];  // 2D array of [x, y] (m/s)
        
        // store the history of body position over time.
        // This is intended to be updated after simulating,
        // although you could give it initial history if you want.
        // intended to be used for plotting.
        this.t_hist = [];  // 1D array of time (s)
        this.p_hist = [];  // 2D array of [x, y] (m)
        this.v_hist = [];  // 2D array of [x, y] (m/s)

    }

    tick(dt) {
        /* tick forward dt seconds time: set the current state of 
        the body to the first element of the expected trajectory,
        and push the current state onto the history. 

        Return true if not able to reach dt, otherwise return false.
        true probably indicates that you should move to the next array
        of bodies in the piecewise array.
        */
        if (this.t.length == 0) {return true;}
        let finalTime = this.t[0] + dt;
        while (this.t[0] < finalTime) {
            this.t_hist.push(this.t.shift())
            this.p_hist.push(this.p.shift())
            this.v_hist.push(this.v.shift())
            if (this.t.length == 0) {return true;}
        }
        return false;
    }
};

let defineSystem = function() {
    /* system is an object representing the current system state.
    The main feature is the pBodies array, which is a piecewise definition
    of the bodies in the system over time.
    The bodies are defined piecewise to allow bodies to be created/destroyed
    when a collision event occurs.
    The arrays of bodies are ordered in the pBodies array increasing in time.

    pBodies[0] ALWAYS contains bodies at the time corresponding to system.t
    That means bodies are shift()ed off the pBodies array when they are used up.
    */
    let system = {
        // t: 0,
        // tol: 1e-8,
        // g: 6.67408e-11,
        pBodies: [[
            new Body('User', 100e2, [-100e9, -100e9], [1, 3.5e4], [0, 0], 100),
            new Body('Sun', 1.98847e30, [0, 0], [0, 0], [0, 0], 695.51e6),
            new Body('Earth', 5.9722e24, [0, 152.10e9], [-29.29e3, 0], [0, 0], 6.371e6),
            new Body('Venus', 4.867e24, [-108.8e9, 0], [0, -35.02e3], [0, 0], 6.0518e6),
            new Body('Ship', 4.867e4, [-108.8e9, 100e9], [0, -200.02e1], [0, 0], 100),
            // new Body('Moon', 7.3477e22, [385e6, 152.10e9], [-29.29e3, 1.022e3], [0, 0], 1.7371e6),
        ]],

        tick: function(dt) {
            /* Move system dt seconds forward in time using the expected values. */
            while (true) {
                let results = new Set([this.pBodies[0].map(body => body.tick(dt))]);
                if (results.has(true)) {this.pBodies.shift();}
                else {break;}
            }
        },

        draw_limit: function() {
            /* given the system, determine the min/max x/y to draw everything. */
            let p = [];
            this.pBodies.forEach( bodies => {
                bodies.forEach( body => {
                    p = p.concat(body.p);
                })
            })
            let x = col_2d(p, 0);
            let y = col_2d(p, 1);
            return [min_1d(x), max_1d(x), min_1d(y), max_1d(y)]
        }
    }
    return system;
}

let deriv_full = function(dydt, y, t, bodies) {
    /*
    q is body position

    qi'' = G*mj / ||qi-qj||**3 * (qi-qj)

    x1 = q
    x2 = q' 

    x1' = x2 = q'
    x2' = q'' = G*mj / ||qi-qj||**3 * (qi-qj)
    */

    let common;
    for (let i = 0; i < bodies.length; i++) {
        // x1' = x2 = q'
        dydt[4*i] = y[4*i+2];  // v in x directino
        dydt[4*i+1] = y[4*i+3];  // v in y direction
        dydt[4*i+2] = 0;  // a in x direction
        dydt[4*i+3] = 0;  // a in y direction
        for (let j = 0; j < bodies.length; j++) {
            if (i == j) { continue; }
            common = (G * bodies[j].m
                      / ((y[4*i] - y[4*j])**2 + (y[4*i+1] - y[4*j+1])**2)**1.5)
            dydt[4*i+2] += common * (y[4*i] - y[4*j]) * -1;
            dydt[4*i+3] += common * (y[4*i+1] - y[4*j+1]) * -1;
        }
    }
}

let get_y0 = function(bodies) {
    /* Return the initial vector for the ODE solver. */ 
    let y0 = [];
    bodies.forEach(function(body){ 
        y0.push(body.p.slice(-1)[0][0]);  // init x
        y0.push(body.p.slice(-1)[0][1]);  // init y
        y0.push(body.v.slice(-1)[0][0]);  // init v in x direction
        y0.push(body.v.slice(-1)[0][1]);  // init v in y direction
    });
    return y0;
}

let collide = function(system) {
    /* If a collision has occured, modify system accordingly and return true.
    
       Modifying system will probably mean creating a new array of bodies
       in the pBodies array where some bodies have been merged.

       Return false if no collision.
    */

   let bodies = system.pBodies.slice(-1)[0];
    let bodyi, bodyj, i_x, i_y, i_r, j_x, j_y, j_r, body_sep;
    for (let i = 0; i < bodies.length; i++) {
        bodyi = bodies[i];
        i_x = bodyi.p.slice(-1)[0][0];
        i_y = bodyi.p.slice(-1)[0][1];
        i_r = bodyi.r;

        for (let j = i+1; j < bodies.length; j++) {

            bodyj = bodies[j];
            j_x = bodyj.p.slice(-1)[0][0];
            j_y = bodyj.p.slice(-1)[0][1];
            j_r = bodyj.r;
            body_sep = ((i_x - j_x)**2  + (i_y - j_y)**2)**0.5;
            if (body_sep <= j_r + i_r) {
                // console.log(bodyi.name, bodyj.name);

                // Create new bodies in system, and merge collided objects.
                
                // debug: reset state.
                if (bodyj.t_hist.length > 0) {
                    bodyj.p[bodyj.p.length] = bodyj.p_hist[0];
                    bodyi.p[bodyi.p.length] = bodyi.p_hist[0];
    
                    bodyj.v[bodyj.v.length] = bodyj.v_hist[0];
                    bodyi.v[bodyi.v.length] = bodyi.v_hist[0];
                }
                else {
                    bodyj.p[bodyj.p.length] = bodyj.p[0];
                    bodyi.p[bodyi.p.length] = bodyi.p[0];
    
                    bodyj.v[bodyj.v.length] = bodyj.v[0];
                    bodyi.v[bodyi.v.length] = bodyi.v[0];
                }
                return true;
            }
        }
    }
    return false;
}

let populate_trajectories = function(system, tIncrease, dt) {
    /* mutate all bodies to update their expected trajectories 
    tIncrease is how many more seconds into the future to simulate.
    Data will be stored in increments of at most dt.
    */

    // Start with the first set of bodies.
    system.pBodies = [system.pBodies[0]];
    // Clear the expected data for the bodies, ready for newly simulated data.
    system.pBodies[0].forEach(body => {
        body.t = [body.t[0]]; 
        body.p = [body.p[0]];
        body.v = [body.v[0]];
    });

    let tSim = system.pBodies[0][0].t[0]; // The time up to which the simulation has been completed.
    let tMax = tSim + tIncrease;  // When simulation reaches tMax, stop.

    while(tSim < tMax) {

        let bodies = system.pBodies.slice(-1)[0];  // Last element is most recent in time.
        let deriv = (dydt, y, t) => deriv_full(dydt, y, t, bodies);
        let y0 = get_y0(bodies); // Initial data
        let integrator = ode45(y0, deriv, 0, dt, {tol: TOL, dtMaxMag: dt}); 
        
        let tInit = tSim;  // Physical time at which ODE starts.

        // Store future trajectories with a time resolution of at least
        // dt, or higher if the ODE solver uses a higher resolution.
        let tEndSub;  // Execute ODE steps until this time is reached. An int multiple of dt.
        let moreStepsRequired;  // True if more steps are required to reach tEndSub.
        let collision;  // True if a collision has occured: integrator needs to be reset.
        while(tSim < tMax) {
            tEndSub = dt*(Math.floor(integrator.t/dt) + 1); // Next highest integer multiple of dt after tSim.
            while(true) {
                moreStepsRequired = integrator.step(tEndSub);
                tSim = integrator.t + tInit;
                bodies.forEach( (body, bodyNum) => {
                    body.t.push(tSim);
                    body.p.push([integrator.y[4*bodyNum], integrator.y[4*bodyNum+1]])
                    body.v.push([integrator.y[4*bodyNum+2], integrator.y[4*bodyNum+3]])
                })
                collision = collide(system);
                if (!moreStepsRequired || collision) {break;}
            }
            if (collision) {break;}  // If collision, break - to force re-init of ODE solver.
        }
    }
    return 1;
}

let plot_orbits = function(system) {
    /* Plot the expected trajectories of the bodies */

    let minx, maxx, miny, maxy;
    [minx, maxx, miny, maxy] = system.draw_limit()
    
    // Apply transformations to make plotting possible in real cartesian coordinates.
    ctx.resetTransform()
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    let xScale = canvas.width/(maxx-minx), yScale = canvas.height/(maxy-miny);
    let scale = Math.min(xScale, yScale); // Same x and y scale
    ctx.transform(1, 0, 0, -1, 0, canvas.height); // Convert to cartesian
    ctx.transform(scale, 0, 0, scale, 0, 0); // scale
    ctx.transform(1, 0, 0, 1, -minx, -miny);  // Shift origin
    ctx.lineWidth = 1/scale;  // Adjust so it's not super thin.

    system.pBodies.forEach( bodies => {
        bodies.forEach(body => {
            drawLine(ctx, body.p);
            ctx.beginPath();
            ctx.arc(body.p[0][0], body.p[0][1], body.r_g, 0, Math.PI*2);
            ctx.fillStyle = "green";
            ctx.fill();
            ctx.closePath();
        });
    })
}

let tick_plot = function(system) {
    let dt = 3600*24;
    populate_trajectories(system, 3600*24*350, dt)
    plot_orbits(system)
    system.tick(dt);
}

let main = function() {
    let system = defineSystem();
    setInterval(tick_plot, 10, system)

}

main()








