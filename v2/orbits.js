// Rotates SVG elements around orbits using animateMotion paths.
 
// theRoot: Should be an SVG element to which the planets 
// will be added.
// thePlanets: Should be all of the SVG elements 
// to be given orbits in the system.
var G= 6.674*Math.pow(10,-11);
var LY2M=1.0*Math.pow(10,16)/1.05702341;
var MSUN=2.0*Math.pow(10,30); //in kg
var MGAL=2.0*Math.pow(10,11);
var RGAL=40000*LY2M;// in m
var SEC2MYEAR=3.16*Math.pow(10,13);
var LYINM40K=3.78*Math.pow(10,20); //40000ly in m

var LYINMONSEC2MYR=8.36*Math.pow(10,-5) //40kly in km x sec in my 
//console.log(Math.pow(10,5));
var V1MSUNINSEC40KLY = G*MSUN/RGAL;//*SEC2MYEAR;
//*Math.pow(10,-13);
//console.log(Math.sqrt(V1MSUNINSEC40KLY*MGAL));


///CLOCK

function createGalaxy(theRoot, theStars, potential) {  
  var galaxy = {
    potential:potential,
    cx: 0,
    cy: 0,
    r: 0, 
    orbitDuration: 1,
    clockwise: true
  };
console.log(galaxy)
  // Private ivars.
  var root = theRoot;
  var stars = theStars;
  var Vinf=105.4;//km/s
  var Rc=0.113341762 //in units of 40KLY//1.39 kpc \\ 4 533.67048 in LY  0.113341762 in 40KLY units  


galaxy.start = function start(d) {
    var system = this;
    for (var i = 0; i < stars.length; ++i) {
      var star = stars[i];
      var path = system.createStarPath(i, star, galaxy.potential);
//      console.log(path,path.getAttribute('id'));
      system.animateGalaxy(star, path.getAttribute('id'), galaxy.potential);
    }
  };

  //  tried to do this without d3 but just could not create the elements 
  // in the right namespaces.
  galaxy.animateGalaxy = function animateGalaxy(star, motionPathId, potential) {
    if (potential == "static"){
      var radinkm=(star.dataset.orbitR/this.cx)*40*LY2M;
      var velkmons=0.000001
    }
    else if (potential == "newtonian"){
	console.log("newtonian");
      var radinkm=(star.dataset.orbitR/this.cx)*40*LY2M;
      var velkmons=Math.sqrt(MGAL*V1MSUNINSEC40KLY/star.dataset.orbitR/this.cx); //radinm;
	if(star.dataset.orbitR/this.cx<0.25){
	    velkmons=velkmons*star.dataset.orbitR/this.cx*4;
	};
      
    }
    else if (potential == "iso")  {
	console.log("iso");
      var radin40KLY=(star.dataset.orbitR/this.cx);//*40*LY2M;
      var RonRc = radin40KLY/Rc
      var velkmons=10*Vinf*Math.sqrt(1.0 - (1.0/RonRc)*Math.atan(RonRc));
    }

    var angvelonmly=  velkmons/(star.dataset.orbitR/this.cx)*LYINMONSEC2MYR;
                    //       console.log("dtheta in MYR",star.dataset.orbitR/this.cx,velkmons,angvelonmly, 2.0*Math.PI/angvelonmly);//*Math.pow(10,-18));
    var animateMotion = d3.select(star).append('animateMotion').attr({
	dur: 2.0*Math.PI/angvelonmly,
	repeatCount: 'indefinite'
    });

//radinm;
//period: 1 sec = 1 My
//1 sec in 
//this.orbitDuration*star.dataset.orbitR/100,
//   console.log(Math.sqrt(G))///*Math.pow(30))))///*Math.pow(10,11)))///(orbit.R*50000/this.cx/LY2M);

    animateMotion.append('mpath').attr('xlink:href', '#' + motionPathId);
    return animateMotion.node();
  };
//      console.log("wtf",animateMotion.getAttribute('dur'))
 
  galaxy.createStarPath = 
  function createStarPath(starIndex, star,potential) {
//      console.log("got here", star.dataset)
//	console.log("cx here",this.cx)
    if (potential == "static"){
	console.log("static");
	var path = d3.select(root).append('path').attr({
	    id: star.id + '-orbit',
	    fill: 'none',
	    stroke: 'black',
	    strokeWidth: 3,
	    'stroke-opacity': 0.7,
	    d: this.makePathDataForOrbit(starIndex, star.dataset.orbitR,star.getAttribute('angle'))
	});
	
	}
      else{
	var path = d3.select(root).append('path').attr({
	    id: star.id + '-orbit',
	    fill: 'none',
	    stroke: 'black',
	    strokeWidth: 3,
	    'stroke-opacity': 0.7,
	    d: this.makePathDataForOrbit(starIndex, star.dataset.orbitR,star.getAttribute('angle'))
	});
      }

    return path.node();
  };
 
  galaxy.makePathDataForOrbit = function makePathDataForOrbit(orbitIndex,
    orbitR,angle) {
    var startPoint = this.positionOnCircle(angle, orbitR);
     this.addDebugCircleD3({x:this.cx,y:this.cy}, "red");
//      console.log(startPoint)
    var rx = orbitR;
    var ry = orbitR;
    var xAxisRotation = 0;
    var largeArcFlag = 1;
    var sweepFlag = this.clockwise ? 1 : 0;
    var vectorToOpposite = {
      x: -2 * (startPoint.x - this.cx),
      y: -2 * (startPoint.y - this.cy),
    };
//    console.log("start point here",this.cx,this.cy,startPoint );
//     this.addDebugCircleD3({
//       x: startPoint.x + vectorToOpposite.x,
//       y: startPoint.y + vectorToOpposite.y
//     },
//    "red");
 
    // Example data path:
    // <path d="M0,400 
    //          a200,200 0 1 0 400,0
    //          a200,200 0 1 0 -400,0"
    //       fill="green" stroke="blue" stroke-width="5" id="redOrbit" />
//     console.log("start point",startPoint.x,startPoint.y,rx,ry);
    var path = 
      // Move command
      'M ' + startPoint.x + ',' + startPoint.y + ' ' + 
    // First arc command
      'a' + rx + ',' + ry + ' ' + 
      xAxisRotation + ' ' + largeArcFlag + ', ' + sweepFlag + ' ' + 
      vectorToOpposite.x + ',' + vectorToOpposite.y + ' ' + 
    // Second arc command
      'a' + rx + ',' + ry + ' ' + 
      xAxisRotation + ' ' + largeArcFlag + ' ,' + sweepFlag + ' ' +
      (-vectorToOpposite.x) + ',' + (-vectorToOpposite.y);
//	  console.log(startPoint.x,startPoint.y)
    return path;
  };

  galaxy.positionOnCircle = function positionOnCircle(angle, radius) {
    var y = radius * Math.sin(+angle);
    var x = radius * Math.cos(+angle);
//      console.log("here",this.cx,this.cy, x, y, x*x+y*y);
    return {x: this.cx + x, y: this.cy + y};
  };
 
  galaxy.addDebugCircleD3 = function addDebugCircleD3(point, color) {
//      console.log("debug",point)
    d3.select(root).append('circle').attr({
      'fill': color,
      'r': 5,
      'cx': point.x,
      'cy': point.y
    });
  };
 
  return galaxy;
}

