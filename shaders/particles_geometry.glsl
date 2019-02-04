layout(triangles) in;
//using this I am at least getting a circle-ish thing
//layout(triangle_strip, max_vertices = 4) out;
//using this I seem to be getting 1/1 ??!!
layout(line_strip, max_vertices = 4) out;

uniform mat4 uv_modelViewProjectionMatrix;
uniform mat4 uv_modelViewInverseMatrix;
uniform float alpha = 0.9;
uniform float uv_fade;
uniform float uv_simulationtimeSeconds;
uniform float drawTrails;

out vec4 color;
out vec2 texcoord;


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

void main()
{

// this will need to be improved to interpolate between values
	vec3 cmap[9]; //Blackbody Colors 0 = 2000 deg, 8 = 10000 deg
    cmap[0]=vec3(1.000, 0.339, 0.174 );
    cmap[1]=vec3(1.000, 0.507, 0.222 );
    cmap[2]=vec3(1.000, 0.720, 0.541  );
    cmap[3]=vec3(1.000, 0.895, 0.809  );
    cmap[4]=vec3(1.000, 0.953, 0.938 );
    cmap[5]=vec3(0.961, 0.953,  1.0 );
    cmap[6]=vec3(0.791, 0.814, 1.0 );
    cmap[7]=vec3(0.540, 0.583, 1.0 );
    cmap[8]=vec3(0.301, 0.359, 1.0 );
	int colIndex = min( max(int(gl_in[0].gl_Position.x/1000. - 2.),0) ,8);
    color = vec4(cmap[colIndex], alpha);
	color = vec4(1.0, 0.0, 0.0, alpha); //red

	//time: this should be the same as simOrbit_geometry.glsl (to synchronize the timing)
	float cosmoTime = gl_in[0].gl_Position.z;
	float Nvals = 4444.;
	float Nsec = 10.;
	float simTime = mod(uv_simulationtimeSeconds,Nsec)*Nvals/Nsec;
	
	//why is this a vec4??
	//vec4 pos = vec4(gl_in[1].gl_Position.x, gl_in[1].gl_Position.y, gl_in[1].gl_Position.z, 1.0);
	vec4 pos = vec4(0., 0., 0., 1.);
    //float rad = 10.0;
	float rad = gl_in[0].gl_Position.y*14950.97871;

	drawSprite(pos, rad, 0);

	//what are the units of radius for drawSprite??
    //if (cosmoTime >= simTime && cosmoTime < (simTime + 1)) {
	//	drawSprite(pos, rad, 0);
	//}
}
