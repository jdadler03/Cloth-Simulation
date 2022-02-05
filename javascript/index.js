function thisIter(S, B){
    return (B+S/B)/2;
}
function babSqrt(S,i,seed){//computes a babylonian square root approximation given S-number and i-iterations
    if(i==0){//seed for computation
    	if(seed==0){
    		var mag = Math.floor(Math.log10(S));//magnitude of number (10^mag)
    		return [S/Math.pow(10, mag+1) + 1.2] * Math.pow(10, mag/2);//rough linear estimation
    	}
    	return seed;
    }
    return thisIter(S, babSqrt(S, i-1, seed));
}

function constraint(P1_,P2_){
	this.P1 = P1_;
	this.P2 = P2_;
}
function vector2(x_,y_){
	this.x = x_;
	this.y = y_;
}
function vMag(V, i=3, seed=0){ return babSqrt(V.x*V.x + V.y*V.y, i, seed);}
function vAdd(V1, V2){ return new vector2(V1.x+V2.x, V1.y+V2.y); }
function vSum(vect, size){
	var sumX=0;
	var sumY=0;
	for(i=0;i<size;i++){
		sumX+=vect[i].x;
		sumY+=vect[i].y;
	}
	return new vector2(sumX, sumY);
}
function vSub(V1,V2){ return new vector2(V1.x-V2.x, V1.y-V2.y); }
function vMply(V1,S){ return new vector2(S*V1.x, S*V1.y)}
function point(X_,isStatic_=0){
	this.X = X_;
	this.isStatic = isStatic_;
	this.A = new vector2(0,0);
}


var ptH = 10;//22;//rows of points quantity
var ptW = 15;//31;//columns of points quantity
var ptD = 30;//distance between points
var ptM = 10;//mass of points
var topMargin = 100;
var leftMargin = 300;
topMargin = topMargin + ptH*ptD;

var grav = new vector2(0, -37034.2);//-20000);//-37034.2);
var pts = new Array;//array of points
var cts = new Array;//array of constraints

for(i=0;i<ptH;i++)//create point array
	for(j=0;j<ptW;j++)
		pts[i*ptW+j] = new point(new vector2(ptD*j , ptD*i));

for (i=651;i<670;i++){
	//pts[i].isStatic=1;
}
for (i=675;i<682;i++){
	//pts[i].isStatic=1;
}
pts[149].isStatic = 1;
pts[135].isStatic = 1;
pts[142].isStatic = 1;
pts[7].isStatic = 1;

var pts0 = pts;

var k = 0;
for(i=0;i<ptH;i++)//generate constraints
	for(j=0;j<ptW;j++){
		if(!(i==ptH-1)) cts[k++] = new constraint(pts[i*ptW+j], pts[i*ptW+j+ptW]);
		if(!(j==ptW-1)) cts[k++] = new constraint(pts[i*ptW+j], pts[i*ptW+j+1]);
	}
var constraintQ = k;//quantity of constraints


function constraintAcc(cts,size){
	var Kconst = 6000;
	for(i=0;i<size;i++){
		var diff = vSub(cts[i].P2.X, cts[i].P1.X);
		var xMag = vMag(diff,2,ptD);
		var xDiff = Math.abs(xMag-ptD);
		if(Math.sign(xMag-ptD)==1) var fMag = xDiff*Kconst;
		else var fMag = 0;
		var F = vMply(diff, fMag/xMag);//force magnitude * normalized difference vector = Force vector
		var A = vMply(F, .99);
		cts[i].P1.A = vAdd(cts[i].P1.A,A);
		cts[i].P2.A = vAdd(cts[i].P2.A, vMply(A,-1));
	}
}

function verlet(dt){
	constraintAcc(cts,constraintQ);
	for(i=0;i<ptH*ptW;i++){
		if(pts[i].isStatic==0) pts[i].X = vAdd(  vSub(vMply(pts[i].X,1.99), vMply(pts0[i].X,.99)) , vMply(pts[i].A, dt*dt)  );
		pts[i].A = new vector2(0,0);
	}
	pts0 = pts;
	constraintAcc(cts,constraintQ);
	for(i=0;i<ptH*ptW;i++){
		if(pts[i].isStatic==0) pts[i].X = vAdd(  vSub(vMply(pts[i].X,1.99), vMply(pts0[i].X,.99)) , vMply(pts[i].A, dt*dt)  );
		pts[i].A = new vector2(0,0);
	}
	pts0 = pts;
	constraintAcc(cts,constraintQ);
	for(i=0;i<ptH*ptW;i++){
		if(pts[i].isStatic==0) pts[i].X = vAdd(  vSub(vMply(pts[i].X,1.99), vMply(pts0[i].X,.99)) , vMply(pts[i].A, dt*dt)  );
		pts[i].A = new vector2(0,0);
	}
	pts0 = pts;
	constraintAcc(cts,constraintQ);
	for(i=0;i<ptH*ptW;i++){
		pts[i].A = vAdd(pts[i].A, grav);
		if(pts[i].isStatic==0) pts[i].X = vAdd(  vSub(vMply(pts[i].X,1.99), vMply(pts0[i].X,.99)) , vMply(pts[i].A, dt*dt)  );
		pts[i].A = new vector2(0,0);
	}
	pts0 = pts;
}

function physics(dt){
	if(dt<.01) verlet(dt);
}


var ctx = document.getElementById("ctx").getContext("2d");
function draw(){
	ctx.clearRect(0, 0, 1000, 800);
	for(i=0;i<ptH;i++)
		for(j=0;j<ptW;j++){
			if(!(i==ptH-1)){
				ctx.beginPath();
				ctx.moveTo(leftMargin + pts[i*ptW+j].X.x, topMargin - pts[i*ptW+j].X.y);
				ctx.lineTo(leftMargin + pts[i*ptW+j+ptW].X.x, topMargin - pts[i*ptW+j+ptW].X.y);
				ctx.stroke();
			}
			if(!(j==ptW-1)){
				ctx.beginPath();
				ctx.moveTo(leftMargin + pts[i*ptW+j].X.x, topMargin - pts[i*ptW+j].X.y);
				ctx.lineTo(leftMargin + pts[i*ptW+j+1].X.x, topMargin - pts[i*ptW+j+1].X.y);
				ctx.stroke();
			}
		}
}

function drawConstraints(){
	for(i=0;i<constraintQ;i++){
		console.log(cts[i]);
		ctx.beginPath();
		ctx.moveTo(leftMargin + cts[i].P1.X.x, topMargin - cts[i].P1.X.y);
		ctx.lineTo(leftMargin + cts[i].P2.X.x, topMargin - cts[i].P2.X.y);
		ctx.stroke();
	}
}


var t0 = null;
function step(t){
	//console.log(ctx.screenX);
	if(!t0) t0 = t;
	physics((t - t0)/1000);
	t0=t;
	draw();
	window.requestAnimationFrame(step);
}
window.requestAnimationFrame(step);

function mouseCoords()
{
  pts[7].X = new vector2(event.clientX - 760,410 - event.clientY);
}



