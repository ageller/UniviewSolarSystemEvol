//NOTE: much of this is overkill because the HZ is never eccentric, 
//but this works and it's a simple build off of the SimOrbit_geometry file
layout(triangles) in;//layout(triangle_strip, max_vertices = 4) out;
layout(triangle_strip, max_vertices = 4) out;
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
out vec4 color;
out vec2 texcoord;
out float texHole;

mat3 rotationMatrix(vec3 axis, float angle) // axis should be normalized
{
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat3(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c);
}
void drawSprite(vec4 position, float radius, float rotation)  // Camera facing square
{  
    vec3 objectSpaceUp = vec3(0, 0, 1);
    vec3 objectSpaceCamera = (uv_modelViewInverseMatrix * vec4(0, 0, 0, 1)).xyz;
    vec3 cameraDirection = normalize(objectSpaceCamera - position.xyz);
    vec3 orthogonalUp = normalize(objectSpaceUp - cameraDirection * dot(cameraDirection, objectSpaceUp));
    vec3 rotatedUp = rotationMatrix(cameraDirection, rotation) * orthogonalUp;
    vec3 side = cross(rotatedUp, cameraDirection);
	rotatedUp=vec3(0,1,0);
	side=vec3(1,0,0);
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

void main()
{
    if(uv_fade<0.001)
        return;

	float TWOPI = 2.0*3.1415926535;
	color = vec4(0,1,0,0.4);
//////////////////////////////////////////////////////////////	

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
		dt = gl_in[2].gl_Position.y/1.e6;
		timeend = simRealtimeend*1000.;
		timestart = simRealtimestart*1000.;
	}
	float usedt = max(dt, simdtmin);
	

//////////////////////////////////////////////////////////////
//Use the data from the three vertices to decode the orbit params
	float lb = mix(gl_in[0].gl_Position.x, gl_in[0].gl_Position.z, (cosmoTime-simTime)/usedt);
	float hb = mix(gl_in[0].gl_Position.y, gl_in[1].gl_Position.x, (cosmoTime-simTime)/usedt);
	texHole=lb/hb;
//////////////////////////////////////////////////////////////
// draw the HZ)
	int drawHZ = 0;
    if ((cosmoTime >= simTime && cosmoTime < (simTime + usedt)) || (simTime >= timeend && cosmoTime >= timeend) || (simTime <= timestart && cosmoTime <= timestart)) {
		drawHZ = 1;
	}
	
	if (drawHZ == 1){
      drawSprite(vec4(0.0,0.0,0.0,1.0),hb*1495.97871,0.0);
	}
	EndPrimitive();
}
