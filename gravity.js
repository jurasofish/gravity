"use strict";

const G = 6.67408e-11; // gravitational constant
const TOL = 1e-8;  // Tolerance for ODE solver.
let canvas = document.getElementById("myCanvas");
let ctx = canvas.getContext("2d");

let windowSize = () => [window.innerHeight, window.innerWidth-200];
[canvas.height, canvas.width] = windowSize();
window.addEventListener('resize', () => {
    [canvas.height, canvas.width] = windowSize();
});

// All coordinates are physical.
let MOUSECLICKED = false; // true if mouse is clicked. probs not perfect.
let MOUSEDOWN = [0, 0];  // Position where mouse was clicked.
let MOUSEPOS = [0, 0];  // Position of mouse hovering.
let MOUSECLICKMATRIX = null;  // transformtion matrix from when mouse was clicked.

let GO = false;  // If true, tick system after every frame.
// GO = true;  // debug
let FINALISEBODIES = false; // True if draft bodies should be finalised in the next tick.
let TICKS = 0;  // How many ticks to move in time while GO is false.

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
    constructor(name, m, p, v, a, r, t=0, draft=false) {
        this.name = name; // Name
        this.m = m; // mass (kg)
        this.a = a; // array of x, y of applied acceleration (m/s/s)
        this.r = r; // Body radius, for for collissions (m)
        // true if the body is still being drawn. This should be set to false
        // on mouseup, indicating that the body has been finalised.
        this.draft = draft;
        
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

    fork() {
        /* Clone the current body and set the clone's initial conditions
           to the current bodie's latest in time expected conditions.
           The cloned conditions are popped from the current body to avoid
           duplication.
           That is, fork the current body.
        */
        let clone = new Body(this.name, 
                             this.m, 
                             this.p.pop(),
                             this.v.pop(), 
                             this.a.slice(), 
                             this.r, 
                             this.t.pop());
        return clone;
    }

    combine(b2) {
        /* return a body which results from the perfectly inelastic collision with body b2 */
        let b1 = this;
        let b1_p = b1.p.slice(-1)[0], b2_p = b2.p.slice(-1)[0];  // body current points.
        let b1_v = b1.v.slice(-1)[0], b2_v = b2.v.slice(-1)[0];  // body current velocities.
        let c_m = b1.m + b2.m;  // Mass of combined bodies.

        // Centroid of the two bodies.
        let c_p = [
            (b1.m * b1_p[0] + b2.m * b2_p[0])/c_m,
            (b1.m * b1_p[1] + b2.m * b2_p[1])/c_m
        ];

        // Conservation of momentum gives final velocity.
        // Assume perfectly inelastic collision.
        let c_v = [
            (b1.m * b1_v[0] + b2.m * b2_v[0])/c_m,
            (b1.m * b1_v[1] + b2.m * b2_v[1])/c_m,
        ];

        let c = new Body(b1.name + ' + ' + b2.name,
                     c_m,
                     c_p,
                     c_v,
                     [0, 0],
                     b1.r + b2.r,  // todo: combine volumes and work out resulting radius assumming constant density
                     b1.t.slice(-1)[0]
        );

        return c;
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

    In any array of bodies, the first body is ALWAYS the user controlled body.
    This will probs change later...
    */
    let system = {
        // t: 0,
        // tol: 1e-8,
        // g: 6.67408e-11,
        pBodies: [[
            new Body('User', 100e2, [-100e9, -100e9], [1, 3.5e2], [0, 0], 100),
            new Body('Sun', 1.98847e30, [0, 0], [0, 0], [0, 0], 695.51e6),
            new Body('Earth', 5.9722e24, [0, 152.10e9], [-29.29e3, 0], [0, 0], 6.371e6),
            new Body('Venus', 4.867e24, [-108.8e9, 0], [0, -35.02e3], [0, 0], 6.0518e6),
            new Body('Ship', 4.867e4, [-108.8e9, 100e9], [0, -200.02e1], [0, 0], 100),
            new Body('Moon', 7.3477e22, [385e6, 152.10e9], [-29.29e3, 1.022e3], [0, 0], 1.7371e6),
        ]],

        tick: function(dt) {
            /* Move system dt seconds forward in time using the expected values. */
            while (true) {
                let results = new Set(this.pBodies[0].map(body => body.tick(dt)));
                if (results.has(true)) {this.pBodies.shift();}
                else {break;}
            }
        },

        untick: function(dt) {
            /* Move system dt seconds back in time using the historical values. */
            // TODO: implement
            console.log('untick not implemented yet');
        },

        draw_limit: function(inputs) {
            /* given the system, determine the min/max x/y to draw everything. 
            The input.follow string is assumed to match one of the names of the bodies,
            in which case a draw limit for that body will be given.
            If no body matches, the entire scene will be used.
            */

            let follow_body = null;

            this.pBodies[0].forEach( body => {
                if (body.name == inputs.follow) {
                    follow_body = body;
                }
            });

            if (follow_body != null) {
                return [
                    follow_body.p[0][0] - 100*follow_body.r, 
                    follow_body.p[0][0] + 100*follow_body.r, 
                    follow_body.p[0][1] - 100*follow_body.r, 
                    follow_body.p[0][1] + 100*follow_body.r, 
                ]
            }

            let p=[];
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
                let bodiesNew = [];  // To be appended to system.pBodies

                // Clone the bodies which aren't involved in the collision.
                bodies.forEach((body, cloneBodyNum) => {
                    // bodiesNew.push(body.fork());  // debug: clone all.
                    if (cloneBodyNum != i && cloneBodyNum != j) {
                        bodiesNew.push(body.fork());  // debug: clone all.
                    }
                });
                bodiesNew.push(bodyi.combine(bodyj));
                system.pBodies.push(bodiesNew)
                return true; // There was a collision - return true.
            }
        }
    }
    return false; // There was no collision - return false.
}

let createDraftBody = function(system, inputs) {
    if (MOUSECLICKED && !GO) {
        let p = MOUSEDOWN;
        let k = inputs.velocity;
        let v = [
            (MOUSEDOWN[0] - MOUSEPOS[0])/k,
            (MOUSEDOWN[1] - MOUSEPOS[1])/k,
        ]
        let t = system.pBodies[0][0].t[0]
        let user = new Body('Userdd', inputs.m, p, v, [0, 0], inputs.r, t, true);
        system.pBodies[0].push(user)
    }
}

let populate_trajectories = function(system, inputs) {
    /* mutate all bodies to update their expected trajectories 
    inputs.lookahead is how many more seconds into the future to simulate.
    Data will be stored in increments of at most inputs.dt.
    */

    // Start with the first set of bodies.
    system.pBodies = [system.pBodies[0]];
    // Clear the expected data for the bodies, ready for newly simulated data.
    system.pBodies[0].forEach(body => {
        body.t = [body.t[0]]; 
        body.p = [body.p[0]];
        body.v = [body.v[0]];
    });

    // Finalise bodies, if inputs triggered it.
    if (FINALISEBODIES) {
        system.pBodies[0].forEach(body => {body.draft=false;});
        FINALISEBODIES=false;
    }

    // Discard the draft bodies.
    system.pBodies[0] = system.pBodies[0].filter(body => !body.draft)

    createDraftBody(system, inputs);  // Create draft body based on current input state.

    let tSim = system.pBodies[0][0].t[0]; // The time up to which the simulation has been completed.
    let tMax = tSim + inputs.lookahead;  // When simulation reaches tMax, stop.

    while(tSim < tMax) {

        let bodies = system.pBodies.slice(-1)[0];  // Last element is most recent in time.
        let deriv = (dydt, y, t) => deriv_full(dydt, y, t, bodies);
        let y0 = get_y0(bodies); // Initial data
        let integrator = ode45(y0, deriv, 0, inputs.dt, {tol: TOL, dtMaxMag: inputs.dt}); 
        
        let tInit = tSim;  // Physical time at which ODE starts.

        // Store future trajectories with a time resolution of at least
        // dt, or higher if the ODE solver uses a higher resolution.
        let tEndSub;  // Execute ODE steps until this time is reached. An int multiple of dt.
        let moreStepsRequired;  // True if more steps are required to reach tEndSub.
        let collision;  // True if a collision has occured: integrator needs to be reset.
        while(tSim < tMax) {
            // Next highest integer multiple of dt after tSim.
            tEndSub = inputs.dt*(Math.floor(integrator.t/inputs.dt) + 1);
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

let plot = function(system, inputs) {
    /* Plot the expected trajectories of the bodies */

    // Clear canvas
    ctx.resetTransform()
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    // Apply transformations to make plotting possible in real cartesian coordinates.
    let minxOrig, maxxOrig, minyOrig, maxyOrig;
    let minx, maxx, miny, maxy;
    [minxOrig, maxxOrig, minyOrig, maxyOrig] = system.draw_limit(inputs)

    // adjust min and max values to zoom in on the centre of the area.
    let x_width = maxxOrig - minxOrig;
    let y_width = maxyOrig - minyOrig;
    let midx = (minxOrig + maxxOrig)/2;
    let midy = (minyOrig + maxyOrig)/2;
    minx = midx - x_width/2/inputs.zoom;
    maxx = midx + x_width/2/inputs.zoom;
    miny = midy - y_width/2/inputs.zoom;
    maxy = midy + y_width/2/inputs.zoom;

    let xScale = canvas.width/(maxx-minx), yScale = canvas.height/(maxy-miny);
    let scale = Math.min(xScale, yScale); // Same x and y scale
    ctx.resetTransform()
    ctx.transform(1, 0, 0, -1, 0, canvas.height); // Convert to cartesian
    ctx.transform(scale, 0, 0, scale, 0, 0); // scale
    ctx.transform(1, 0, 0, 1, -minx, -miny);  // Shift origin
    ctx.lineWidth = 1/scale;  // Adjust so it's not super thin.

    system.pBodies.forEach( (bodies, bodiesIdx) => {
        bodies.forEach(body => {
            if (bodiesIdx == 0) {  // Only draw bodies for the first set of bodies.
                let r_g = Math.log(body.r) * inputs.bodysize;  // graphical radius.
                ctx.beginPath();
                ctx.arc(body.p[0][0], body.p[0][1], r_g, 0, Math.PI*2);
                ctx.fillStyle = "green";
                ctx.fill();
                ctx.closePath();
            }
            drawLine(ctx, body.p);  // Draw paths for all sets of bodies.
        });
    });

    if (MOUSECLICKED) {
        // Plot 
        drawLine(ctx, [MOUSEDOWN, MOUSEPOS]);
    };
}

let get_inputs = function() {
    /* return inputs and a boolean flag indicating whether they are valid,
       and set error message */

    let errFlag = false;

    let set_error = (m) => document.getElementById("error").innerHTML = m;
    set_error('');  // Clear it

    let inputs = {
        m: Number(document.getElementById("mass-text").value),
        r: Number(document.getElementById("radius-text").value),
        // The ODE solver does NOT like high precision float inputs.
        dt: Math.round(Number(document.getElementById("dt-text").value)*3600*24),
        lookahead: Math.round(Number(document.getElementById("lookahead-text").value)*3600*24),
        
        bodysize: Number(document.getElementById("bodysize-text").value),
        zoom: Number(document.getElementById("zoom-text").value),
        velocity: Number(document.getElementById("velocity-text").value),
        follow: document.getElementById("follow-text").value,
    }

    Object.keys(inputs).forEach((key) => {
        if (typeof inputs[key] === "number" && isNaN(inputs[key])) {
            set_error(key + ' is invalid');
            errFlag = true;
            return;
        }
    });
    if (errFlag) {return [inputs, errFlag];};
    
    // example to update the slider and text box programmatically.
    // document.getElementById("mass-text").value = (inputs.m * 1.01).toExponential(2);
    // document.getElementById("mass-text").oninput();

    return [inputs, errFlag];

}

let tick_plot = function(system) {
    let inputs, errFlag;
    [inputs, errFlag] = get_inputs()
    if (errFlag) {return;}

    // Based on current mouse click and drag status,
    // modify the bodies to include the user body.

    populate_trajectories(system, inputs)
    plot(system, inputs)
    if (GO || TICKS > 0) {
        system.tick(inputs.dt);
        TICKS = Math.max(0, TICKS-1);
    }
    if (TICKS < 0) {
        system.untick(inputs.dt);
        TICKS = Math.min(0, TICKS+1);
    }
}

let main = function() {
    let system = defineSystem();
    setInterval(tick_plot, 10, system)
}

// Return current matrix to transform from canvas to physical coords
let getTransMatrix = () => ctx.getTransform().invertSelf();

let transformCoords = function(p, mat=getTransMatrix()) {
    return [
        p[0] * mat.a + p[1] * mat.c + mat.e,
        p[0] * mat.b + p[1] * mat.d + mat.f
    ];
}

canvas.addEventListener('mousedown', e => {
    let x = e.clientX - canvas.offsetLeft;
    let y = e.clientY - canvas.offsetTop;
    MOUSECLICKMATRIX = getTransMatrix();
    GO = false;
    [x, y] = transformCoords([x, y], MOUSECLICKMATRIX);
    MOUSEDOWN = [x, y];
    MOUSEPOS = [x, y];
    MOUSECLICKED = true;
});

canvas.addEventListener('mouseup', e => {
    if (MOUSECLICKED) {
        // GO = true;
    }
    MOUSECLICKED = false;
    FINALISEBODIES = true;
});

canvas.addEventListener('mouseout', e => {
    MOUSECLICKED = false;
});

canvas.addEventListener('mousemove', e => {
    let x = e.clientX - canvas.offsetLeft;
    let y = e.clientY - canvas.offsetTop;
    if (MOUSECLICKMATRIX == null) {
        MOUSECLICKMATRIX = getTransMatrix();
    }
    [x, y] = transformCoords([x, y], MOUSECLICKMATRIX);
    MOUSEPOS = [x, y];
});

document.addEventListener("keydown", event => {
    console.log(event.keyCode);
    if(event.target.tagName.toLowerCase() != 'body') {
        return;
    };

    switch(event.keyCode) {
        case 32:  // space bar
            GO = !GO; break;
        case 37: // left arrow
            TICKS--; break;
        case 39: // right arrow
            TICKS++; break;
        case 67: // c key
            MOUSECLICKED = false; break;
        default:
            break;
    }
});

main()
