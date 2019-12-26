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

let draw_limit = function(bodies) {
    /* given the bodies, determine the min/max x/y to draw everything. */
    let p = [];
    bodies.forEach( body => {
        p = p.concat([body.p], body.p_exp);
    })
    //bodies.map(body => p = p.concat([body.p], body.p_exp));
    let x = col_2d(p, 0);
    let y = col_2d(p, 1);
    return [min_1d(x), max_1d(x), min_1d(y), max_1d(y)]
}

let Body = class{
    /* A body is a physical object which produces and is affected by gravity. */
    constructor(name, m, p, v, a, r, t=0, t_hist=[], p_hist=[], 
                v_hist=[], t_exp=[], p_exp=[], v_exp=[]) {
        this.name = name; // Name
        this.m = m; // mass (kg)
        this.p = p; // array of x, y of position (m)
        this.v = v; // array of x, y of velocity (m/s)
        this.a = a; // array of x, y of applied acceleration (m/s/s)
        this.r = r; // Body radius, for for collissions (m)
        this.r_g = Math.log(r) * 5e8; // Body radius, for graphics.
        this.t = t; // Time at which the body is at point p.
        
        // store the history of body position over time.
        // This is intended to be updated after simulating,
        // although you could give it initial history if you want.
        // intended to be used for plotting.
        this.t_hist = t_hist;  // 1D array of time (s)
        this.p_hist = p_hist;  // 2D array of [x, y] (m)
        this.v_hist = v_hist;  // 2D array of [x, y] (m/s)
        
        // Store the expected future position over time.
        // Similar to t_hist and p_hist.
        this.t_exp = t_exp;  // 1D array of time (s)
        this.p_exp = p_exp;  // 2D array of [x, y] (m)
        this.v_exp = v_exp;  // 2D array of [x, y] (m/s)

    }

    tick() {
        /* tick one step forward in time: set the current state of 
        the body to the first element of the expected trajectory,
        and push the current state onto the history. 
        */
        this.t_hist.push(this.t)
        this.p_hist.push(this.p)
        this.v_hist.push(this.v)

        this.t = this.t_exp.shift()
        this.p = this.p_exp.shift()
        this.v = this.v_exp.shift()
    }
};

let defineBodies = function() {

    let bodies = [
        new Body('User', 100e2, [-100e9, -100e9], [1, 3.5e4], [0, 0], 100),
        new Body('Sun', 1.98847e30, [0, 0], [0, 0], [0, 0], 695.51e6),
        new Body('Earth', 5.9722e24, [0, 152.10e9], [-29.29e3, 0], [0, 0], 6.371e6),
        new Body('Venus', 4.867e24, [-108.8e9, 0], [0, -35.02e3], [0, 0], 6.0518e6),
        new Body('Ship', 4.867e4, [-108.8e9, 100e9], [0, -200.02e2], [0, 0], 100),
        // new Body('Moon', 7.3477e22, [385e6, 152.10e9], [-29.29e3, 1.022e3], [0, 0], 1.7371e6),
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

let update_bodies = function(bodies, nBodies, t_solved, y_solved) {
    for (let bodyNum = 0; bodyNum < nBodies; bodyNum++) {
        let body = bodies[bodyNum]
        body.t_exp = t_solved.map(x => x + body.t);
        body.p_exp = y_solved.map(x => [x[4*bodyNum], x[4*bodyNum+1]]);
        body.v_exp = y_solved.map(x => [x[4*bodyNum+2], x[4*bodyNum+3]]);
    }
    /*
    for (let i = 0; i < nBodies; i++) {
        for (let j = 0; j < nBodies; j++) {
            let body_sep = ((bodies[i].p[0] - bodies[j].p[0])**2 
                           + (bodies[i].p[1] - bodies[j].p[1])**2)**0.5;
            console.log(bodies[i].name + ' ' + bodies[j].name + ' ' + body_sep)
        }
    }
    */
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

    // Now update the bodies.
    update_bodies(bodies, nBodies, t_solved, y_solved)
}

let plot_orbits = function(bodies) {
    /* Plot the expected trajectories of the bodies */

    let data = []
    bodies.forEach(body => {
        var trace = {
            x: col_2d(body.p_exp, 0),
            y: col_2d(body.p_exp, 1),
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
    // Plotly.newPlot('plot', data, layout);

    let minx, maxx, miny, maxy;
    [minx, maxx, miny, maxy] = draw_limit(bodies)
  
    ctx.resetTransform()
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    let xScale = canvas.width/(maxx-minx), yScale = canvas.height/(maxy-miny);
    let scale = Math.min(xScale, yScale); // Same x and y scale
    ctx.transform(1, 0, 0, -1, 0, canvas.height); // Convert to cartesian
    ctx.transform(scale, 0, 0, scale, 0, 0); // scale
    ctx.transform(1, 0, 0, 1, -minx, -miny);  // Shift origin
    ctx.lineWidth = 1/scale;  // Adjust so it's not super thin.

    bodies.forEach(body => {
        drawLine(ctx, body.p_exp);
        ctx.beginPath();
        ctx.arc(body.p[0], body.p[1], body.r_g, 0, Math.PI*2);
        ctx.fillStyle = "green";
        ctx.fill();
        ctx.closePath();
    });
}

let tick_plot = function(bodies, nBodies) {
    populate_trajectories(bodies, nBodies, 3600*24*350, 3600*24)
    plot_orbits(bodies)
    bodies.forEach(body => body.tick())
}

let main = function() {
    let bodies, nBodies;
    [bodies, nBodies] = defineBodies()
    setInterval(tick_plot, 10, bodies, nBodies)

}

main()








