/*
 * LICENCE
 * copyright 2014 ~ ****
 * Some rights reserved.
 * Author: HUFANGYUAN
 * Released under CC BY-NC
*/
#version 330 core

flat in int fType;
in vec4 fPos;
in float f3;
in float f4;
out vec4 color;
uniform mat4 fMvpInv;
uniform int flag; 
uniform float range;

const vec4 white = vec4(1.f, 1.f, 1.f, 1.f);
const vec4 red = vec4(1.f, 0.f, 0.f, 1.f);
const vec4 green = vec4(0.f, 1.f, 0.f, 1.f);
const vec4 blue = vec4(0.f, 0.f, 1.f, 1.f);

void showCircle();
void outPres();
void outVort();

void main()
{
	//showCircle();

	switch (flag) {
	case 1: //pres
		outPres();
		break;
	case 2: //vorticity
		outVort();
		break;
	default:
		break;
	}

	return;
}

void showCircle() {
	vec2 pixPos = gl_PointCoord.st * vec2(2.f, -2.f) - vec2(1.f, -1.f);
	if (length(pixPos) >= 1.f) discard;
	if (abs(pixPos.x - pixPos.y)>0.2 && abs(pixPos.x + pixPos.y)>0.2) discard;
}

void outPres() {
	float tmp = f4 / range;
	if (fType == 0 || fType == 1) {
		//color = vec4(0., length(fVel)/5, 0., 1.f);
		if (tmp >= 0.f && tmp < 0.5f)
			color = 2.*((0.5f - tmp)* blue + tmp* green);
		else if (tmp >= 0.5f && tmp < 1.f)
			color = 2.*((1.f - tmp)* green + (tmp - 0.5f)* red);
		else if (tmp < 0.f)
			color = blue;
		else
			color = red;
	}
	/*
	else if(fType == 1)
	{
	color = vec4(1.f, 0.f, 0.f, 1.f);
	}
	*/
	else {
		color = vec4(0.f, 0.f, 1.f, 1.f);
	}
	color.a = 1.f;
}

void outVort() {
	float tmp = f3 / range;
	if(fType == 0 || fType == 1) {
		if (tmp >= -1.f && tmp < 0.f)
			color = white - ( (0.f - tmp)* (red + green) );
		else if (tmp >= 0.f && tmp < 1.f)
			color = white - ( (tmp - 0.f)* (green + blue) );

		else if (tmp < -1.f)
			color = blue;
		else
			color = red;
	}
	else {
		color = vec4(0.f, 0.f, 1.f, 1.f);
	}
	color.a = 1.f;
}