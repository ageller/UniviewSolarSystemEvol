in vec2 texcoord;
in vec4 color;
in float texHole;

out vec4 fragColor;

void main()
{
    float texpos= length(2*(texcoord-vec2(0.5,0.5)));
	if (texpos>1.0||texpos<texHole){
	  discard;
	}
	
    fragColor = color;
}
