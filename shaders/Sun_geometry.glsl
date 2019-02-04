layout(triangles) in;
layout(triangle_strip, max_vertices = 128) out;
//layout(line_strip, max_vertices = 128) out;
uniform float divLon = 6.; //max_vertices
uniform float divLat = 6.; //max_vertices
uniform float alpha = 1.0; //for Sun
uniform mat4 uv_projectionMatrix;
uniform mat4 uv_modelViewMatrix;
uniform mat4 uv_modelViewInverseMatrix;
uniform mat4 uv_modelViewProjectionMatrix;
uniform mat4 uv_normalMatrix;
uniform int uv_simulationtimeDays;
uniform float uv_simulationtimeSeconds;
uniform float uv_fade;
out vec2 texcoord;
out vec4 color;
//
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

uniform float uLevel = 2;
uniform float uLayers = 3;


// axis should be normalized
mat3 rotationMatrix(vec3 axis, float angle)
{
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat3(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c);
}

void drawSprite(vec4 position, float radius, float rotation)
{
    vec3 objectSpaceUp = vec3(0, 0, 1);
    vec3 objectSpaceCamera = (uv_modelViewInverseMatrix * vec4(0, 0, 0, 1)).xyz;
    vec3 cameraDirection = normalize(objectSpaceCamera - position.xyz);
    vec3 orthogonalUp = normalize(objectSpaceUp - cameraDirection * dot(cameraDirection, objectSpaceUp));
    vec3 rotatedUp = rotationMatrix(cameraDirection, rotation) * orthogonalUp;
    vec3 side = cross(rotatedUp, cameraDirection);
    texcoord = vec2(0, 1);
	gl_Position = uv_modelViewProjectionMatrix * vec4(position.xyz + radius * (-side + rotatedUp), 1);
	EmitVertex();
    texcoord = vec2(0, 0);
	gl_Position = uv_modelViewProjectionMatrix * vec4(position.xyz + radius * (-side - rotatedUp), 1);
	EmitVertex();
    texcoord = vec2(1, 1);
	gl_Position = uv_modelViewProjectionMatrix * vec4(position.xyz + radius * (side + rotatedUp), 1);
	EmitVertex();
    texcoord = vec2(1, 0);
	gl_Position = uv_modelViewProjectionMatrix * vec4(position.xyz + radius * (side - rotatedUp), 1);
	EmitVertex();
	EndPrimitive();
}

void ProduceVertex( float s, float t, float rad, vec3 V0, vec3 V01, vec3 V02) { 
	vec3 v = V0 + s*V01 + t*V02; 	
	v = normalize(v); 
	vec4 n = vec4(v, 1.0); 
	vec4 tnorm = normalize( uv_normalMatrix * n );  // the transformed normal
	vec4 ECposition = uv_modelViewMatrix * vec4( (rad*v), 1. ); 
	gl_Position = uv_projectionMatrix * ECposition; 
	EmitVertex( );
}

float linterp(vec2 p1, vec2 p2, float xout){
	float m = (p2[1] - p1[1])/(p2[0] - p1[0]);
	float b = p2[1] - m*p2[0];
	
	return m*xout + b;

}


void main()
{
    if(uv_fade<0.001)
        return;

	float PI = 3.1415926535;

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
		dt = gl_in[2].gl_Position.y/1.e6;
		timeend = simRealtimeend*1000.;
		timestart = simRealtimestart*1000.;
	}
	float usedt = max(dt, simdtmin);
	
//	//time: this should be the same as simOrbit_geometry.glsl (to synchronize the timing)
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
// colors
	vec3 cmap[9]; //Blackbody Colors 0 = 2000 deg, 8 = 10000 deg
//original
//    cmap[0]=vec3(1.000, 0.339, 0.174 );
//    cmap[1]=vec3(1.000, 0.507, 0.222 );
//    cmap[2]=vec3(1.000, 0.720, 0.541 );
//    cmap[3]=vec3(1.000, 0.895, 0.809 );
//    cmap[4]=vec3(1.000, 0.953, 0.938 );
//    cmap[5]=vec3(0.961, 0.953,  1.0 );
//    cmap[6]=vec3(0.791, 0.814, 1.0 );
//    cmap[7]=vec3(0.540, 0.583, 1.0 );
//    cmap[8]=vec3(0.301, 0.359, 1.0 );
//exagerated
	cmap[0]=vec3(1.000, 0.139, 0.074 ); //2000
    cmap[1]=vec3(1.000, 0.207, 0.074 ); //3000
    cmap[2]=vec3(1.000, 0.520, 0.441 ); //4000
    cmap[3]=vec3(1.000, 0.695, 0.509 ); //5000
    cmap[4]=vec3(1.000, 0.953, 0.938 ); //6000
    cmap[5]=vec3(0.861, 0.953,  1.0 ); //7000
    cmap[6]=vec3(0.691, 0.714, 1.0 ); //8000
    cmap[7]=vec3(0.440, 0.583, 1.0 ); //9000
    cmap[8]=vec3(0.201, 0.359, 1.0 ); //10000
//	int colIndex = min( max(int(gl_in[0].gl_Position.x/1000. - 2.),0) ,8);
//	color = vec4(cmap[colIndex], alpha);
//	color = vec4(1.0, 0.0, 0.0, alpha); //red
//instead I want to linearly interpolate along this cmap (so that the colors don't jump from one to the other)
	float teff = mix(gl_in[0].gl_Position.x, gl_in[0].gl_Position.z, (cosmoTime-simTime)/usedt);
	float cval = teff/1000. - 2.;
	int I1 = int(min( max(floor(cval),0) ,8));
	int I2 = int(min( max(ceil(cval),0) ,8));
	if (I1 == I2){
		color = vec4(cmap[I1], alpha);
	} else {
		float red   = linterp(vec2(I1, cmap[I1].x), vec2(I2, cmap[I2].x), cval);
		float green = linterp(vec2(I1, cmap[I1].y), vec2(I2, cmap[I2].y), cval);
		float blue  = linterp(vec2(I1, cmap[I1].z), vec2(I2, cmap[I2].z), cval);
		color = vec4(red, green, blue, alpha);
	}
//	int colIndex = min( max(int(cval),0) ,8);
//	color = vec4(cmap[colIndex], alpha);


	
	vec4 pos;

//////////////////////////////////////////////////////////////
// how should we draw the Sun?  The Sprite is the best option.	
	int doIcosasphere = 0;
	int doSprite = 1;
	int doAMGspere = 0;
	
	//for icosasphere (getting errors unless these are defined, even if doIcosaspher == 0)
	vec3 V0s[8];
	V0s[0] = vec3(0., 0., 1.);
	V0s[1] = vec3(1., 0., 0.);
	V0s[2] = vec3(0., 0., -1.);
	V0s[3] = vec3(-1., 0., 0.);
	V0s[4] = vec3(0., 0., 1.);
	V0s[5] = vec3(1., 0., 0.);
	V0s[6] = vec3(0., 0., -1.);
	V0s[7] = vec3(-1., 0., 0.);

	vec3 V1s[8];
	V1s[0] = vec3(1., 0., 0.);
	V1s[1] = vec3(0., 0., -1.);
	V1s[2] = vec3(-1., 0., 0.);
	V1s[3] = vec3(0., 0., 1.);
	V1s[4] = vec3(1., 0., 0.);
	V1s[5] = vec3(0., 0., -1.);
	V1s[6] = vec3(-1., 0., 0.);
	V1s[7] = vec3(0., 0., 1.);

	vec3 V2s[8];
	V2s[0] = vec3(0., 1., 0.);
	V2s[1] = vec3(0., 1., 0.);
	V2s[2] = vec3(0., 1., 0.);
	V2s[3] = vec3(0., 1., 0.);
	V2s[4] = vec3(0., -1., 0.);
	V2s[5] = vec3(0., -1., 0.);
	V2s[6] = vec3(0., -1., 0.);
	V2s[7] = vec3(0., -1., 0.);

//////////////////////////////////////////////////////////////
// draw the sun (?)
	int drawOrb = 0;
//    if ((simTime >= cosmoTime && simTime < (cosmoTime + dt)) || (simTime == timeend && cosmoTime >= timeend) || (simTime == timestart && cosmoTime <= timestart)) {
    if ((cosmoTime >= simTime && cosmoTime < (simTime + usedt)) || (simTime >= timeend && cosmoTime >= timeend) || (simTime <= timestart && cosmoTime <= timestart)) {

		drawOrb = 1;
	}
	
	if (drawOrb == 1){
		//draw sphere from triangles

//		float rad = gl_in[0].gl_Position.y*1495.97871;
		float rad = mix(gl_in[0].gl_Position.y, gl_in[1].gl_Position.x, (cosmoTime-simTime)/usedt)*1495.97871;
	
		if (doSprite == 1){
			pos = vec4(0., 0., 0., 1.);
			drawSprite(pos, rad, 0);
		}
		
//following http://web.engr.oregonstate.edu/~mjb/cs519/Handouts/geometry_shaders.2pp.pdf
		if (doIcosasphere == 1) {
			for (int i=0; i<8; i++){
			
				vec3 V01 = V1s[i] - V0s[i];
				vec3 V02 = V2s[i] - V0s[i];
				vec3 V0 = V0s[i];
				
				float dt = 1. / float( uLayers );
				float t_top = 1.;
				for( int it = 0; it < uLayers; it++ ) { 	
					float t_bot = t_top - dt; 
					float smax_top = 1. - t_top; 
					float smax_bot = 1. - t_bot;
				
					int nums = it + 1; 
					float ds_top = smax_top / float( nums - 1 ); 
					float ds_bot = smax_bot / float( nums );
				
					float s_top = 0.; 
					float s_bot = 0.;
			
					for( int is = 0; is < nums; is++ ) { 
						ProduceVertex( s_bot, t_bot, rad, V0, V01, V02);
						ProduceVertex( s_top, t_top, rad, V0, V01, V02);
				
						s_top += ds_top; 
						s_bot += ds_bot; 
					}
					ProduceVertex( s_bot, t_bot, rad, V0, V01, V02);
					EndPrimitive();
			
					t_top = t_bot; 
					t_bot -= dt;
				}
			}
		}

	

	
//following http://stackoverflow.com/questions/28776292/trying-to-make-a-sphere-with-triangle-strips-on-opengl
		if (doAMGspere == 1){
			float x;
			float y;
			float z;
		
			for (int i=0;i<divLon;i++) {
				for (int j=0;j<=divLat;j++){
			
					if (mod(i,2) == 0){
						if (abs(z) != rad || (j != 0 && j != divLat)){
							x = sin(PI*j/divLat)*cos(2.*PI*i/divLon)*rad;	
							y = sin(PI*j/divLat)*sin(2.*PI*i/divLon)*rad;	
							z = cos(PI*j/divLat)*rad;
							pos = vec4(x, y, z, 1.0);
							gl_Position = uv_modelViewProjectionMatrix * pos;
							EmitVertex();					
						}
						if (abs(z) != rad && j != 0) {
							x = sin(PI*j/divLat)*cos(2.*PI*(i+1)/divLon)*rad;	
							y = sin(PI*j/divLat)*sin(2.*PI*(i+1)/divLon)*rad;	
							z = cos(PI*j/divLat)*rad;
							pos = vec4(x, y, z, 1.0);
							gl_Position = uv_modelViewProjectionMatrix * pos;
							EmitVertex();							
						}
					}
					if (mod(i,2) != 0){
						if (abs(z) != rad || (j != 0 && j != divLat)){
							x = sin(PI*(divLat - j)/divLat)*cos(2.*PI*i/divLon)*rad;	
							y = sin(PI*(divLat - j)/divLat)*sin(2.*PI*i/divLon)*rad;	
							z = cos(PI*(divLat - j)/divLat)*rad;
							pos = vec4(x, y, z, 1.0);
							gl_Position = uv_modelViewProjectionMatrix * pos;
							EmitVertex();					
						}
						if (abs(z) != rad && j != 0) {
							x = sin(PI*(divLat - j)/divLat)*cos(2.*PI*(i+1)/divLon)*rad;	
							y = sin(PI*(divLat - j)/divLat)*sin(2.*PI*(i+1)/divLon)*rad;	
							z = cos(PI*(divLat - j)/divLat)*rad;
							pos = vec4(x, y, z, 1.0);
							gl_Position = uv_modelViewProjectionMatrix * pos;
							EmitVertex();							
						}					
					}
				}
			}
			
			EndPrimitive();

		}		
	}



}
