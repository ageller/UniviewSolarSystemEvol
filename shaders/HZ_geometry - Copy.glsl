//NOTE: much of this is overkill because the HZ is never eccentric, 
//but this works and it's a simple build off of the SimOrbit_geometry file
layout(triangles) in;
//layout(triangle_strip, max_vertices = 4) out;
layout(triangle_strip, max_vertices = 128) out;
//layout(line_strip, max_vertices = 128) out;
uniform float Nvertices = 64; //max_vertices
uniform mat4 uv_projectionMatrix;
uniform mat4 uv_modelViewMatrix;
uniform mat4 uv_modelViewInverseMatrix;
uniform mat4 uv_modelViewProjectionMatrix;
uniform int uv_simulationtimeDays;
uniform float uv_simulationtimeSeconds;
uniform float uv_fade;
//
uniform int doBack;
//uniform float simLoop;
uniform float simRealtime;
//uniform float simTimefac;
uniform float simRealtimestart;
uniform float simRealtimeend;
//uniform float simRealtimeepoch;
uniform float simMasstimestart;
uniform float simMasstimeend;
//uniform float simShowstartonly;
uniform float simdtmin;
//out vec2 texcoord;
out vec4 color;

void main()
{
    if(uv_fade<0.001)
        return;

	float TWOPI = 2.0*3.1415926535;
	
//////////////////////////////////////////////////////////////	
//colors: probably some better way to do this, but for now let's do the simple thing
	float cval = (gl_in[2].gl_Position.z);
	if (cval == 1){ //Mercury
		color = vec4(1.0, 1.0, 1.0, 1.0); //white
	}
	if (cval == 2){ //Venus
		color = vec4(0.7, 0.5, 0.6, 1.0); //brown
	}
	if (cval == 3){ //Earth
		color = vec4(0.0, 0.0, 1.0, 1.0); //blue
	}
	if (cval == 4){ //Mars
		color = vec4(1.0, 0.0, 0.0, 1.0); //red
	}
	if (cval == 5){ //Jupiter
		color = vec4(1.0, 0.5, 0.0, 1.0); //orange
	}
	if (cval == 6){ //Saturn
		color = vec4(0.5, 0.5, 0.5, 1.0); //gray
	}
	if (cval == 7){ //Uranus
		color = vec4(0.0, 1.0, 1.0, 1.0); //aqua
	}
	if (cval == 8){ //Neptune
		color = vec4(1.0, 1.0, 0.0, 1.0); //yellow
	}
	if (cval == 9){ //Pluto
		color = vec4(0.7, 0.2, 0.9, 1.0); //purple
	}
	if (cval >= 10){ //Habitabe Zone
		color = vec4(0.0, 0.0, 1.0, 0.7); //blue
	}
	vec3 globalUp = vec3(0,1,0);

//////////////////////////////////////////////////////////////
//define the time (sim will run from year 0 to year 13000, each year representing one Myr)
	float dayfract = uv_simulationtimeSeconds/(24.0*3600.0);//0.5*2.0*3.14*(time)/(sqrt(a.x*a.x*a.x/3347937656.835192));
    float years_0 = (uv_simulationtimeDays + dayfract)/365.2425 + 1970.;
    float cosmoTime = clamp(years_0,0.0,13800.0);     
	float simTime = gl_in[1].gl_Position.z;
	float dt = 1.;
	float timeend = simMasstimeend;
	float timestart = simMasstimestart;
	if (simRealtime == 1 ){ 
		simTime = gl_in[2].gl_Position.x/1.e6;
//		dt = gl_in[2].gl_Position.y/1.e6;
		dt = gl_in[2].gl_Position.y/1.e6;
		timeend = simRealtimeend*1000.;
		timestart = simRealtimestart*1000.;
	}
	float usedt = max(dt, simdtmin);
	
//	//time steps (full sim takes Nsec and the number of frames in Nvals)
//	//set this up first assuming equal time steps in mass loss
//	float cosmoTime = gl_in[1].gl_Position.z;
//	float dt = 1.;
//	float dayfract = uv_simulationtimeSeconds/(24.0*3600.0);
//	float epoch = simRealtimeepoch;
//	int days = uv_simulationtimeDays;
//	float currentyr = (days - epoch + dayfract)*365.2425;
//	float timestart = simMasstimestart;
//	float timeend = simMasstimeend;
//	float Nfac = 1.e6/simTimefac;
//	if (simLoop == 0){
//		Nfac = 0.01*simTimefac;//to get roughly the same time scale
//	}
//	//if real time is desired
//	if (simRealtime == 1 ){ 
//		cosmoTime = gl_in[2].gl_Position.x;
//		dt = gl_in[2].gl_Position.y;
//		timestart = simRealtimestart*1.e9;
//		timeend = simRealtimeend*1.e9;
//		Nfac = 1.e6/simTimefac;
//		if (simLoop == 0){
//			Nfac = 200.*simTimefac;//to get roughly the same time scale
//		}
//	}
//	float Nvals = timeend - timestart;
//	float simTime = mod(currentyr, Nfac)*Nvals/Nfac + timestart;
//	//if looping is desabled
//	if (simLoop == 0){
//		simTime = currentyr * Nfac + timestart; //yr * Nfac
//		if (simTime > timeend){
//			simTime = timeend;
//		}
//		if (simTime < timestart){
//			simTime = timestart;
//		}
//	}

//////////////////////////////////////////////////////////////
//Use the data from the three vertices to decode the orbit params
//	float a = gl_in[0].gl_Position.x; //AU
//	float a2 = gl_in[0].gl_Position.y; //AU
	float a = mix(gl_in[0].gl_Position.x, gl_in[0].gl_Position.z, (cosmoTime-simTime)/usedt);
	float a2 = mix(gl_in[0].gl_Position.y, gl_in[1].gl_Position.x, (cosmoTime-simTime)/usedt);
	float e = 0.0;
	float i = 0.0;
	float OMEGA = 0.0;
	float omega = 0.0;
	
//////////////////////////////////////////////////////////////
// draw the orbits (?)
	int drawOrb = 0;
//    if ((simTime >= cosmoTime && simTime < (cosmoTime + dt)) || (simTime == timeend && cosmoTime >= timeend) || (simTime == timestart && cosmoTime <= timestart)) {
    if ((cosmoTime >= simTime && cosmoTime < (simTime + usedt)) || (simTime >= timeend && cosmoTime >= timeend) || (simTime <= timestart && cosmoTime <= timestart)) {
		drawOrb = 1;
	}
	
	if (drawOrb == 1){
		//Transform

		vec3 b = vec3(-i,omega,OMEGA);
		vec3 c = cos(b);
		vec3 s = sin(b);
		vec4 pos; 
		vec4 pos1;
		vec4 pos2;
		vec4 pos3;	
		a = a* 1495.97871;
		a2 = a2* 1495.97871;
		float E=0.0;
		vec3 P;
		P.x = c.z*c.y - s.z*c.x*s.y;
		P.y = c.z*s.y + s.z*c.x*c.y;
		P.z = s.z*s.x;
		P=vec3(-P.x,-P.y,-P.z);
		vec3 Q;
		Q.x = -s.z*c.y - c.z*c.x*s.y;
		Q.y = -s.z*s.y + c.z*c.x*c.y;
		Q.z = s.x*c.z;
		Q = vec3(-Q.x,-Q.y,-Q.z);
		
        //max_vertices from above minus 1
		float Ntheta = Nvertices-1.;
		float dTheta = TWOPI/Ntheta;
		
		for (int i=0;i<=Ntheta;i++) {
//		for (int i=int(Ntheta);i>=0;i--) {
			int j = i;
			if (doBack == 1) {
				j = int(Ntheta) - i;
			}
			E = j*dTheta;
			pos1 = vec4(a * (cos(E) - e) * P + a * sqrt(1.0 - e * e) * sin(E) * Q, 1.0);
			gl_Position = uv_modelViewProjectionMatrix * pos1;
			EmitVertex();
				
				
			E = j*dTheta;
			pos1 = vec4(a2 * (cos(E) - e) * P + a2 * sqrt(1.0 - e * e) * sin(E) * Q, 1.0);
			gl_Position = uv_modelViewProjectionMatrix * pos1;
			EmitVertex();
				
				
		}
	}
	EndPrimitive();
}
