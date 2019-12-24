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
        let x = point[0];
        let y = point[1];
        ctx.lineTo(x,y);
    });
    ctx.stroke();
    ctx.closePath();
}

let Body = class{
    /* A body is a physical object which produces and is affected by gravity. */
    constructor(name, m, p, v, a, 
                t_hist=nj.float64([]), p_hist=nj.float64([]), 
                t_exp=nj.float64([]), p_exp=nj.float64([])) {
        this.name = name; // Name
        this.m = m; // mass (kg)
        this.p = p; // array of x, y of position (m)
        this.v = v; // array of x, y of velocity (m/s)
        this.a = a; // array of x, y of applied acceleration (m/s/s)
        
        // store the history of body position over time.
        // This is intended to be updated after simulating,
        // although you could give it initial history if you want.
        // intended to be used for plotting.
        this.t_hist = t_hist;  // 1D float64 ndarray of time (s)
        this.p_hist = p_hist;  // 2D float64 ndarray of [x, y] (m)
        
        // Store the expected future position over time.
        // Similar to t_hist and p_hist.
        this.t_exp = t_exp;  // 1D float64 ndarray of time (s)
        this.p_exp = p_exp;  // 2D float64 ndarray of [x, y] (m)

    }
};

let defineBodies = function() {

    let bodies = [
        new Body('Sun', 1.98847e30, [0, 0], [0, 0], [0, 0]),
        new Body('Earth', 5.9722e24, [0, 152.10e9], [-29.29e3, 0], [0, 0]),
        new Body('Venus', 4.867e24, [-108.8e9, 0], [0, -35.02e3], [0, 0]),
    ]

    let nBodies = bodies.length  // Number of bodies

    return [bodies, nBodies];
    
}

let deriv_full = function(dydt, y, t, bodies, nBodies) {

    /*

    q is body position

    qi'' = G*mj / ||qi-qj||**3 * (qi-qj)

    x1 = q
    x2 = q' 

    x1' = x2 = q'
    x2' = q'' = G*mj / ||qi-qj||**3 * (qi-qj)

    */
    
    for (let i = 0; i < nBodies; i++) {

        // x1' = x2 = q'
        dydt[4*i] = y[4*i+2];  // v in x directino
        dydt[4*i+1] = y[4*i+3];  // v in y direction
        
        dydt[4*i+2] = 0;  // a in x direction
        dydt[4*i+3] = 0;  // a in y direction

        for (let j = 0; j < nBodies; j++) {
            if (i == j) { continue; }
            // console.log(i, j)

            let body_sep = ((y[4*i] - y[4*j])**2 
                           + (y[4*i+1] - y[4*j+1])**2)**0.5;
        
            dydt[4*i+2] += G * bodies[j].m / body_sep**3 * (y[4*i] - y[4*j]) * -1;
            dydt[4*i+3] += G * bodies[j].m / body_sep**3 * (y[4*i+1] - y[4*j+1]) * -1;

        }
    }
}

let get_y0 = function(bodies) {
    /* Return the initial vector for the ODE solver. */ 
    let y0 = []; // Initial data
    bodies.forEach(function(body){ 
        y0.push(body.p[0]);  // init x
        y0.push(body.p[1]);  // init y
        y0.push(body.v[0]);  // init v in x direction
        y0.push(body.v[1]);  // init v in y direction
    });
    return y0;
}

let populate_trajectories = function(bodies, nBodies, tmax, dt) {
    /* mutate all bodies to update their expected trajectories */

    let y0 = get_y0(bodies); // Initial data
    let deriv = function(dydt, y, t) {
        return deriv_full(dydt, y, t, bodies, nBodies)
    };
    let integrator = ode45(y0, deriv, 0, dt, {tol: TOL});
    
    // Integrate up to tmax in steps of exactly dt.
    let y_solved = [], t_solved = [];
    while(integrator.t < tmax) {
        integrator.steps(Infinity, integrator.t + dt)
        t_solved.push(integrator.t);
        y_solved.push(integrator.y.slice());
    }

    t_solved = nj.float64(t_solved)
    y_solved = nj.float64(y_solved)

    // Now update the bodies.
    for (let bodyNum = 0; bodyNum < nBodies; bodyNum++) {
        bodies[bodyNum].t_exp.push = t_solved;
        bodies[bodyNum].p_exp = nj.stack(
            [y_solved.pick(null, 4*bodyNum), y_solved.pick(null, 4*bodyNum+1)],
            -1
        );
    }
}

let plot_orbits = function(bodies) {
    /* Plot the expected trajectories of the bodies */

    let data = []
    bodies.forEach(body => {
        var trace = {
            x: body.p_exp.pick(null, 0).tolist(),
            y: body.p_exp.pick(null, 1).tolist(),
            mode: 'lines+markers',
            type: 'scatter'
        };
        data.push(trace)
    });
    var layout = {
        yaxis: {
        scaleanchor: "x",
        },
    }
    Plotly.newPlot('plot', data, layout);
}

let main = function() {
    let bodies, nBodies;
    [bodies, nBodies] = defineBodies()
    populate_trajectories(bodies, nBodies, 3600*24*350, 3600*24*7*4)
    plot_orbits(bodies)
}

main()

let xy_pts_test = [
    [100, 100],
    [150, 150],
    [200, 300],
    [400, 230]
];
drawLine(ctx, xy_pts_test)







