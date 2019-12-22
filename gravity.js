function drawLine(ctx, pts) {
    ctx.beginPath();
    pts.forEach(point => {
        var x = point[0];
        var y = point[1];
        ctx.lineTo(x,y);
    });
    ctx.stroke();
    ctx.closePath();
}

G = 6.67408e-11 // gravitational constant
var canvas = document.getElementById("myCanvas");
var ctx = canvas.getContext("2d");
height = 400
width = 600
canvas.height = height
canvas.width = width

var xy_pts_test = [
    [100, 100],
    [150, 150],
    [200, 300],
    [400, 230]
];
drawLine(ctx, xy_pts_test)

// x and y coordinates of bodies (m)
xy = [
    [0, 0],  // sun, defined at origin.
    [0, 152.10e9],  // Earth
    [-108.8e9, 0] // Venus
]

// x and y component of velocities of bodies (m/s)
v = [
    [0, 0],  // Sun
    [-29.29e3, 0],  // Earth
    [0, -35.02e3] // Venus
]

// masses of bodies in (kg)
m = [
    1.98847e30,  // Sun
    5.9722e24,  // Earth
    4.867e24  // Venus
]

nBodies = xy.length  // Number of bodies

var y0 = [] // Initial data
for (let i = 0; i < nBodies; i++) {
    y0.push(xy[i][0])  // init x
    y0.push(xy[i][1])  // init y
    y0.push(v[i][0])  // init v in x direction
    y0.push(v[i][1])  // init v in y direction
}
console.log(y0)

/*

q is body position

qi'' = G*mj / ||qi-qj||**3 * (qi-qj)

x1 = q
x2 = q' 

x1' = x2 = q'
x2' = q'' = G*mj / ||qi-qj||**3 * (qi-qj)

*/

var deriv = function(dydt, y, t) {
    
    for (let i = 0; i < nBodies; i++) {

        // x1' = x2 = q'
        dydt[4*i] = y[4*i+2]  // v in x directino
        dydt[4*i+1] = y[4*i+3]  // v in y direction
        
        dydt[4*i+2] = 0  // a in x direction
        dydt[4*i+3] = 0  // a in y direction

        for (let j = 0; j < nBodies; j++) {
            if (i == j) { continue; }
            // console.log(i, j)

            body_sep = ((y[4*i] - y[4*j])**2 
                        + (y[4*i+1] - y[4*j+1])**2)**0.5
        
            dydt[4*i+2] += G * m[j] / body_sep**3 * (y[4*i] - y[4*j]) * -1
            dydt[4*i+3] += G * m[j] / body_sep**3 * (y[4*i+1] - y[4*j+1]) * -1
            
            bnlah = 1;

        }
    }
    // console.log("dydt")
    // console.log(dydt)
}


// Initialize: 
// var integrator = ode45(y0, deriv, 0, 60, {"tol": 1e-18})
var integrator = ode45(y0, deriv, 0, 3600)
 
// Integrate up to tmax: 
var n = 0, tmax = 3600*24*365.256363004, t = [], y_out = []
while( integrator.step( tmax ) ) {
  // Store the solution at this timestep: 
  t.push( integrator.t )
  y_out.push( integrator.y.slice() )
  n++

  console.log("t; y")
  // console.log(t)
  // console.log(y_out)
}
console.log(n)
console.log(t)
console.log(y_out)

var data = []
for (let i = 0; i < nBodies; i++) {
    var x_pts = [], y_pts = []
    for (let l = 0; l < t.length; l++) {
        x_pts.push(y_out[l][4*i+0])
        y_pts.push(y_out[l][4*i+1])
    }
    var trace = {
        x: x_pts,
        y: y_pts,
        mode: 'lines+markers',
        type: 'scatter'
      };
    data.push(trace)
}
var layout = {
    yaxis: {
      scaleanchor: "x",
    },
}
Plotly.newPlot('plot', data, layout);