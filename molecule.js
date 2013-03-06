var atoms = [];
var bonds = [];
var xCoords = [];
var yCoords = [];
var zCoords = [];
var size = [];
var colors = [];
var elements = [];
var bondCount = 0;		
var connections = [];
var lineStack = [];
var zStack = [];
var atomDepth = 0;
var rad = Math.PI/180;		
var Xang, Yang, Zang;
var xAng = Xang = 0;
var yAng = Yang = 0;
var xAng = Zang = 0;
var doRotate = false;
var rotateByAxis = false;
var prevX;
var prevY;

var colorTable = [];
colorTable['C'] = "333333";
colorTable['H'] = "CCCCCC";
colorTable['O'] = "CC0033";
colorTable['N'] = "3399FF";
colorTable['HE'] = "FFCCFF";
colorTable['LI'] = "996600";
colorTable['B'] = "009900";
colorTable['F'] = "CCCC00";
colorTable['NA'] = "006699";
colorTable['MG'] = "336600";
colorTable['AL'] = "000000";
colorTable['SI'] = "CCCC00";
colorTable['P'] = "FF6600";
colorTable['S'] = "FFCC00";
colorTable['CL'] = "009900";
colorTable['CA'] = "666666";
colorTable['TI'] = "666666";
colorTable['CR'] = "666666";
colorTable['MN'] = "666666";
colorTable['FE'] = "FF6600";
colorTable['NI'] = "996600";
colorTable['CU'] = "996600";
colorTable['ZN'] = "996600";
colorTable['BR'] = "996600";
colorTable['AG'] = "666666";
colorTable['I'] = "3399FF";
colorTable['BA'] = "FF6600";
colorTable['AU'] = "CCCC00";
colorTable['unknown'] = "FF66FF";
			
var radiusTable = [];
radiusTable['C'] = 67;
radiusTable['H'] = 53;
radiusTable['O'] = 48;
radiusTable['N'] = 56;
radiusTable['P'] = 98;
radiusTable['F'] = 42;
radiusTable['S'] = 88;
radiusTable['BR'] = 94;
radiusTable['CL'] = 79;
radiusTable['I'] = 115;

var mol = new Object();
mol.atoms = [];
mol.bonds = [];
mol.connections = [];

var canvas;
var ctx;

function trace(str) {
     console.log(str);
}
			
function sorter(a, b) {
	return a[0]*1000000 - b[0]*1000000;
}

function spinAtoms() {
	zStack = [];
	renderAtoms();
	rotate(1, 1, 1);
}

function rotate(xRot, yRot, zRot) {
	rotateByAxis = true;
	zStack = [];
	Xang += xRot * rad;
	Yang += yRot * rad;
	Zang += zRot * rad;
			
	for (var i = 0; i < this.atoms.length; i++) {
		var atom = mol.atoms[i];
		// rotation
		var cosy = Math.cos(Yang);
		var cosx = Math.cos(Xang);
		var siny = Math.sin(Yang);
		var sinx = Math.sin(Xang);
		var sinz = Math.sin(Zang);
		var cosz = Math.cos(Zang);
		var Xpos1 = atom.coords[0];
		var Ypos1 = atom.coords[1] * cosx + atom.coords[2] * sinx;
		var Zpos1 = atom.coords[1] * -sinx + atom.coords[2] * cosx;
		var Xpos2 = Xpos1 * cosy - Zpos1* siny;
		var Ypos2 = Ypos1;
		var Zpos2 = Xpos1 * siny + Zpos1 * cosy;
		var Xpos3 = Xpos2 * cosz - Ypos2 * sinz;
		var Ypos3 = Ypos2 * cosz + Xpos2 * sinz;
		var Zpos3 = Zpos2;
		var depth = 1/(1-(Zpos3/250));
		atom.y = 2 * Ypos3 * depth + 200;
		atom.x = 2 * Xpos3 * depth + 200;
		// scale depending on depth
		atom.scale = (0.5*Zpos3+depth+atom.size)*0.6/100;
		atom.radius = atom.scale*25;
		atom.Zpos = Zpos3;
		atom.rad = (atom.Zpos*depth+atom.size) * .15;
		// add atom to zStack depending on its depth
		zStack.push([depth, atom]);
	}
	renderAtoms();
}

function addAtom (atomX, atomY, atomZ, size, element) {
	var thisAtom = new Object();
	thisAtom.neighbors = [];
	thisAtom.lines = [];
	thisAtom.coords = [];
	thisAtom.coords[0] = atomX;
	thisAtom.coords[1] = atomY;
	thisAtom.coords[2] = atomZ;
	thisAtom.size = size;
	thisAtom.element = element;
	thisAtom.name = "a" + atomDepth;
	thisAtom.id = atomDepth;
	mol.atoms.push(thisAtom);
	atomDepth++;
}

function addLines(lines) {
	for (var i = 0; i < lines.length; i++) {
		var thisBond = new Object();
		var id = mol.bonds.length;
		thisBond.name = "b" + id;
		thisBond.id = id;
		mol.bonds.push(thisBond);
		atomDepth++;
		var atom1 = mol.atoms[lines[i][0]];
		var atom2 = mol.atoms[lines[i][1]];
		var x1 = atom1.coords[0];
		var y1 = atom1.coords[1];
		var z1 = atom1.coords[2];
		var x2 = atom2.coords[0];
		var y2 = atom2.coords[1];
		var z2 = atom2.coords[2];
		var mag = Math.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
		mol.connections.push([atom1, atom2, mag]);
		atom1.lines.push(thisBond);
		atom2.lines.push(thisBond);
		atom1.neighbors.push(atom2);
		atom2.neighbors.push(atom1);
	}
	checkBonds();
}

function checkBonds () {
	// loop through all atoms to see if each is "happy"
	for (var i = 0; i < mol.atoms.length; i++) {
		var atom = mol.atoms[i];
	}
	// loop through all bonds
	var bondStrength;
	for (i = 0; i < mol.connections.length; i++) {
		// if we find a bond that binds 2 unhappy atoms, it needs to be double
		var a1 = mol.connections[i][0];
		var a2 = mol.connections[i][1];
		bondStrength = 1;
		a1.happy = isAtomHappy(a1);
		a2.happy = isAtomHappy(a2);
		if (!a1.happy && !a2.happy) {
			a1.neighbors.push(a2);
			a2.neighbors.push(a1);
			bondStrength++;
		}
		if (!a1.happy && !a2.happy) {
			a1.neighbors.push(a2);
			a2.neighbors.push(a1);
			bondStrength++;
		}
	}
}

function getMouse(event, canvasElem){ 
	var mouse = {x:0, y:0}; 
	var isTouch = false; 

	// touch screen events 
	if (event.touches){ 
		isTouch = true; 
		if (event.touches.length){ 
			mouse.x = parseInt(event.touches[0].pageX); 
			mouse.y = parseInt(event.touches[0].pageY); 
		} 
	}else{ 
		// mouse events 
		mouse.x = parseInt(event.clientX); 
		mouse.y = parseInt(event.clientY); 
	} 

	// accounts for border 
	mouse.x -= canvasElem.clientLeft; 
	mouse.y -= canvasElem.clientLeft; 

	// parent offsets 
	var par = canvasElem; 
	while (par !== null) { 
		if (isTouch){ 
			// touch events offset scrolling with pageX/Y 
			// so scroll offset not needed for them 
			mouse.x -= parseInt(par.offsetLeft); 
			mouse.y -= parseInt(par.offsetTop); 
		}else{ 
			mouse.x += parseInt(par.scrollLeft - par.offsetLeft); 
			mouse.y += parseInt(par.scrollTop - par.offsetTop); 
		} 

		par = par.offsetParent; 
	} 

	return mouse; 
} 

function isAtomHappy(a) {
	if (a.element == "C" && a.neighbors.length < 4) {
		isHappy = false;
	} else if (a.element == "N" && a.neighbors.length < 3) {
		isHappy = false;
	} else if ((a.element == "O" || a.element == "S") && a.neighbors.length < 2) {
		isHappy = false;
	} else {
		isHappy = true;
	}
	return isHappy;
}

function renderAtoms() {
	if (!rotateByAxis) {
		for (var i = 0; i < mol.atoms.length; i++) {
			var atom = mol.atoms[i];
			//
			var Xang = xAng * rad;
			var Yang = yAng * rad;
			// rotation
			var cosy = Math.cos(Yang);
			var cosx = Math.cos(Xang);
			var siny = Math.sin(Yang);
			var sinx = Math.sin(Xang);
			var tempZ = atom.coords[2] * cosy - atom.coords[0] * siny;
			var Xpos = atom.coords[2] * siny + atom.coords[0] * cosy;
			var Zpos = atom.coords[1] * sinx + tempZ * cosx;
			var Ypos = atom.coords[1] * cosx - tempZ * sinx;
			var depth = 1/(1-(Zpos/250));
			atom.x = 2 * Xpos * depth + 200;
			atom.y = 2 * Ypos * depth + 200;
						
			atom.Zpos = Zpos;
			atom.rad = (atom.Zpos*depth+atom.size) * .15;
			// scale depending on depth
			atom.scale = (0.5*Zpos+depth+atom.size)*0.6/100;
			// add atom to zStack depending on its depth
			zStack[i] = [depth, atom];
		}
	}
	zStack.sort(sorter);
	lineStack = [];
	var tempStack = [];
	// rearrange zStack
	for (i = 0; i < zStack.length; i++) {
		atom = zStack[i][1];
		tempStack.push(["atom", atom]);
		for (var j = 0; j < atom.lines.length; j++) {
			var bond = atom.lines[j];
			if (lineStack[bond.id] == undefined) {
				tempStack.push(["bond", bond]);
				lineStack[bond.id] = true;
			}
		}
	}
			
	zStack = null;
	zStack = tempStack;
	drawMolecule();
}

function drawMolecule() {
	ctx.clearRect(0, 0, 500, 500);
	for (var i = 0; i < zStack.length; i++) {
		if (zStack[i][0] == "atom") {
			var atom = zStack[i][1];
			//trace(zStack[i][1].x + ": " + atom.scale + ": " + atom.size);
			var grdOffset = 3;
      		var grd=ctx.createRadialGradient(atom.x - grdOffset,atom.y - grdOffset,grdOffset,atom.x,atom.y,25*atom.scale);
			grd.addColorStop(0,  "white");
			grd.addColorStop(0.9,colorTable[atom.element]);
			grd.addColorStop(1, 'rgba(102,102,102,0)');

			// Fill with gradient
			ctx.fillStyle=grd;
			//ctx.fillRect(0, 0, 500, 350);
			ctx.arc(atom.x, atom.y, atom.radius, 0, Math.PI*2, true);
			ctx.fill();
			ctx.closePath();
						
		} else if (zStack[i][0] == "bond") {
			var bond = zStack[i][1];
			var atom1 = mol.connections[bond.id][0];
			var atom2 = mol.connections[bond.id][1];

    		var v0, v1, x1, y1;
			if (atom1.Zpos < atom2.Zpos) {
				v0 = atom1;
				v1 = atom2;
				x1 = atom2.x;
				y1 = atom2.y;
			} else {
				v0 = atom2;
				v1 = atom1;
				x1 = atom1.x;
				y1 = atom1.y;
			}
			// make sure bonds walk around smoothly on perimeter
			var zDiff = v1.Zpos - v0.Zpos;
			var percent = 1 - zDiff/mol.connections[bond.id][2];
			var radius = atom1.radius * percent;
			var xSide = v0.x - v1.x;
			var ySide = v0.y - v1.y;
			var hyp = Math.sqrt((xSide*xSide) + (ySide*ySide));
			var theta = Math.abs(Math.asin(xSide/hyp));
			//trace("theta: " + theta);
			var xDir;
			var yDir;
			v0.x < v1.x ? xDir = 1 : xDir = -1;
			v0.y < v1.y ? yDir = 1 : yDir = -1;
			// hack to prevent lines from not showing up
			if (theta == 0) {
				theta = 0.001;
			}
			// end hack
		
			var xPoint = v0.x + (radius*Math.sin(theta))*xDir;
			var yPoint = v0.y + (radius*Math.cos(theta))*yDir;
			// put bond in place and scale it
			var x = xPoint;
			var y = yPoint;
			var scaleX = (x1 - xPoint)/100;
			var scaleY = -(yPoint - y1)/100;

			ctx.beginPath();
    		ctx.moveTo(x, y);
    		ctx.lineTo(x + 100*scaleX, y + 100*scaleY);
    		ctx.closePath();
    		ctx.stroke();

		}
	}
}

function init() {

	var xmlhttp = new XMLHttpRequest();
	xmlhttp.addEventListener("load", onXMLLoad, false);
	xmlhttp.open("GET", "caffeine.cml", false);
	xmlhttp.send("");
	var xmlDoc = xmlhttp.responseXML;
	canvas = document.getElementById('myCanvas');
    ctx = canvas.getContext('2d');
	spinAtoms();

	canvas.addEventListener('mousemove', molMove, false);
	canvas.addEventListener('mouseout', molOut);
	canvas.addEventListener('mousedown', molDown);
	canvas.addEventListener('mouseup', molUp);
}


function molOut(ev) {

	doRotate=false;

}

function molDown(ev) {

	prevX = getMouse(ev, canvas).x;
	prevY = getMouse(ev, canvas).y;
	doRotate = true;
}

function molUp(ev) {

	doRotate = false;

}

function molMove(ev) {

	var locX = prevX;
	var locY = prevY;
	var xDiff;
	var yDiff;

	if(doRotate){

		prevX = getMouse(ev, canvas).x;
		prevY = getMouse(ev, canvas).y;
		xDiff = prevX - locX;
		yDiff = prevY - locY;
		yAng -= xDiff;
		xAng -= yDiff;

		rotate(yDiff, -xDiff, 0);
	}
}

function parseXML(xmlDoc) {
	// get atoms from XML
	atoms = [];
	for (i = 0; i < xmlDoc.getElementsByTagName("atom").length; i++) {
		node = xmlDoc.getElementsByTagName("atom")[i];
		the_id = node.getAttribute("id");
		trace(the_id);
		atom = new Object();
		atoms.push(atom);
		atom.id = Number(the_id);
		for (ii = 0; ii < node.getElementsByTagName("string").length; ii++) {
			var newnode = node.getElementsByTagName("string")[ii];
			if (newnode.getAttribute("title")) {
				atom.element = newnode.firstChild.nodeValue;
			}
		} 
		coordNode = node.getElementsByTagName("coordinate3")[0];
		atom.coords = [];
		atom.coords = coordNode.firstChild.nodeValue.split(" ");
		trace(atom.coords);
	}
	// get bonds from XML
	bonds = [];
	bondCount = 0;
	for (i = 0; i < xmlDoc.getElementsByTagName("list").length; i++) {
		node = xmlDoc.getElementsByTagName("list")[i];
		if (node.getAttribute("title") == "connections") {
			for (ii = 0; ii < node.getElementsByTagName("list").length; ii++) {
				var newnode = node.getElementsByTagName("list")[ii];
				first_atom = Number(newnode.getAttribute("id") - 1);
				trace(first_atom);
				connections = [];
				connections = newnode.firstChild.nodeValue.split(" ");
				trace(connections);
				for (b = 0; b < connections.length; b++) {
				// if the connection is to an atom whose ID isn't 0
					if (connections[b] > 0) {
						// add the connection to a 2-dimensional array [atom1 ID,atom2 ID]
						bonds.push([first_atom, connections[b] - 1]);
						bondCount++;
					}
				}
			} 
		}
					
	}
	getReady();
}

function onXMLLoad(event) {
	parser = new DOMParser();
	xmlDoc = parser.parseFromString(event.target.responseText,"text/xml");
	parseXML(xmlDoc);
}

function getReady() {
	xCoords = [];
	yCoords = [];
	zCoords = [];
	size = [];
	colors = [];
	lines = [];
	elements = [];
	// loop through atom array from XML obj.
	for (var i = 0; i < atoms.length; i++) {
		var atom = atoms[i];
		// read its coordinates into x, y, and z arrayCoords
		xCoords.push(Number(atom.coords[0]));
		yCoords.push(Number(atom.coords[1]));
		zCoords.push(Number(atom.coords[2]));
		// determine color and size depending on element
		try {
			size.push(200);
			//size.push(radiusTable[atom.element.toUpperCase()]*3.5);
		} catch(e) {
			size.push(200);
		}
		var the_color;
		try {
			the_color = colorTable[atom.element.toUpperCase()];
		}catch(e) {
			the_color = colorTable['unknown'];
		}
		colors.push(the_color);
		elements.push(atom.element.toUpperCase());
	}
	// loop through bond array from XML obj.;
	for (i = 0; i < bondCount; i++) {
		var bondExists = false;
		// since CML file lists each bond twice (once for each atom it binds),
		// test if the bond already exists
		for (var j = 0; j < lines.length; j++) {
			if (bonds[i][0] == lines[j][1] && bonds[i][1] == lines[j][0])
			{
				bondExists = true;
			}
		}
		// if it doesn't, add it to the lines array
		if (!bondExists)
		{
			lines.push(bonds[i]);
		}
	}
	var x_sum = 0;
	var y_sum = 0;
	var z_sum = 0;
	for (i = 0; i < xCoords.length; i++)
	{
		yCoords[i] *=  -1;
		x_sum +=  Number(xCoords[i]);
		y_sum +=  Number(yCoords[i]);
		z_sum +=  Number(zCoords[i]);
	}
	// get average x,y, and z
	var x_avg = x_sum / xCoords.length;
	var y_avg = y_sum / yCoords.length;
	var z_avg = z_sum / zCoords.length;
	adjustCoordinates(x_avg,xCoords);
	adjustCoordinates(y_avg,yCoords);
	adjustCoordinates(z_avg,zCoords);
	// loop through x array and add each atom
	atomDepth = 0;
	for (i = 0; i < xCoords.length; i++)
	{
		addAtom(xCoords[i] * 18,yCoords[i] * 18,zCoords[i] * 18,size[i] / 2,elements[i]);
	}
	addLines(lines);
			
};

function adjustCoordinates(avg,arr)
{
	if (Math.abs(avg) > 0)
	{
		for (var i = 0; i < arr.length; i++)
		{
			arr[i] -=  avg;
		}
	}
}			