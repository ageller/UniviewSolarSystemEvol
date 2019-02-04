//NOTE: this is essentially a copy of simOrbits__geometry.glsl to only show the initial values.  
// I'm sure that I could be more clever about this and not require two nearly identical files!
layout(triangles) in;
//layout(triangle_strip, max_vertices = 4) out;
//layout(triangle_strip, max_vertices = 128) out;
layout(line_strip, max_vertices = 128) out;
uniform float Nvertices = 128.; //max_vertices
uniform mat4 uv_projectionMatrix;
uniform mat4 uv_modelViewMatrix;
uniform mat4 uv_modelViewInverseMatrix;
uniform mat4 uv_modelViewProjectionMatrix;
uniform int uv_simulationtimeDays;
uniform float uv_simulationtimeSeconds;
uniform float uv_fade;
//out vec2 texcoord;
out vec4 color;


void main()
{
    if(uv_fade<0.001)
        return;

	float TWOPI = 2.0*3.1415926535;
	
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
	if (cval >= 10){ //Habitable Zone
		color = vec4(0.0, 0.0, 1.0, 0.7); //blue
	}
	vec3 globalUp = vec3(0,1,0);
	
	//Use the data from the three vertices to decode the orbit params
	float a = gl_in[0].gl_Position.x; //AU
	float a2 = gl_in[2].gl_Position.x; //AU
	float e = gl_in[0].gl_Position.y;
	float i = gl_in[1].gl_Position.x*TWOPI/360.;
//	float omega = (gl_in[0].gl_Position.z- gl_in[1].gl_Position.y )*TWOPI/360.;
	float omega = gl_in[1].gl_Position.y*TWOPI/360.; // argument of periapsis
	float OMEGA = gl_in[0].gl_Position.z*TWOPI/360. ; // Longitude of ascending nodes
 //	float omega = (gl_in[0].gl_Position.z + gl_in[1].gl_Position.y )*TWOPI/360.;

	//time steps (full sim takes Nsec and the number of frames in Nvals)
	float cosmoTime = gl_in[1].gl_Position.z;
	float Nvals = 4444.;
	float Nsec = 10.;
	float simTime = mod(uv_simulationtimeSeconds,Nsec)*Nvals/Nsec;
//    if (cosmoTime >= simTime && cosmoTime < (simTime + 1)) {
    if (cosmoTime == 0) {
		//Transform

		vec3 b = vec3(-i,OMEGA, omega);
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
		
		if (a2 <= 0){
			for (int i=0;i<=Ntheta;i++) {
				E = i*dTheta;
				pos = vec4(a * (cos(E) - e) * P + a * sqrt(1.0 - e * e) * sin(E) * Q, 1.0);
				gl_Position = uv_modelViewProjectionMatrix * pos;
				EmitVertex();
			}
		}
	}
	EndPrimitive();
}
