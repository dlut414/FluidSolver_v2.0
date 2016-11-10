/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
#define GLSL(version, shader)  "#version " #version "\n" #shader

const GLchar* fragment0 = GLSL(330 core,

uniform mat4 fMvpInv;
uniform float sRangeMax;
uniform float sRangeMin;

flat in int fType;
in vec4 fPos;
in float fS1;
in float fS2;
out vec4 color;

const vec4 white = vec4(1.f, 1.f, 1.f, 1.f);
const vec4 red = vec4(1.f, 0.f, 0.f, 1.f);
const vec4 green = vec4(0.f, 1.f, 0.f, 1.f);
const vec4 blue = vec4(0.f, 0.f, 1.f, 1.f);

void showCircle();
void paintRGB();
void paintRB();

void main() {
	paintRGB();
	return;
}

void showCircle() {
	vec2 pixPos = gl_PointCoord.st * vec2(2.0f, -2.0f) - vec2(1.0f, -1.0f);
	if (length(pixPos) >= 1.0f) discard;
	if (abs(pixPos.x - pixPos.y)>0.2f && abs(pixPos.x + pixPos.y)>0.2f) discard;
}

void paintRGB() {
	float range = sRangeMax - sRangeMin;
	float s_normalized = (fS1 - sRangeMin) / range;
	if (fType == 2) {
		color = vec4(0.1f, 0.1f, 0.1f, 1.0f);
	}
	else if (fType == 3) {
		color = vec4(1.0f, 1.0f, 0.0f, 1.0f);
	}
	else if (fType == 4) {
		color = vec4(0.0f, 1.0f, 1.0f, 1.0f);
	}
	else if(fType == 0) {
		if (s_normalized >= 0.0f && s_normalized < 0.5f)
			color = 2.0f * ((0.5f - s_normalized)* blue + s_normalized* green);
		else if (s_normalized >= 0.5f && s_normalized < 1.0f)
			color = 2.0f * ((1.0f - s_normalized)* green + (s_normalized - 0.5f)* red);
		else if (s_normalized < 0.0f)
			color = blue;
		else
			color = red;
	}
	color.a = 1.0f;
}

void paintRB() {
	float range = sRangeMax - sRangeMin;
	float s_normalized = 2.0f* (fS1 - sRangeMin) / range - 1.0f;
	if (fType == 2) {
		color = vec4(0.1f, 0.1f, 0.1f, 1.0f);
	}
	else if (fType == 3) {
		color = vec4(1.0f, 1.0f, 0.0f, 1.0f);
	}
	else if (fType == 4) {
		color = vec4(0.0f, 1.0f, 1.0f, 1.0f);
	}
	else if (fType == 0) {
		if (s_normalized >= -1.0f && s_normalized < 0.0f)
			color = white - ((0.0f - s_normalized)* (red + green));
		else if (s_normalized >= 0.0f && s_normalized < 1.0f)
			color = white - ((s_normalized - 0.0f)* (green + blue));
		else if (s_normalized < -1.0f)
			color = blue;
		else
			color = red;
	}
	color.a = 1.0f;
}

);

const GLchar* fragment1 = GLSL(330 core,

uniform mat4 fMvpInv;
uniform float sRangeMax;
uniform float sRangeMin;

flat in int fType;
in vec4 fPos;
in float fS1;
in float fS2;
out vec4 color;

const vec4 white = vec4(1.f, 1.f, 1.f, 1.f);
const vec4 red = vec4(1.f, 0.f, 0.f, 1.f);
const vec4 green = vec4(0.f, 1.f, 0.f, 1.f);
const vec4 blue = vec4(0.f, 0.f, 1.f, 1.f);

void showCircle();
void paintRGB();
void paintRB();

void main() {
	paintRGB();
	return;
}

void showCircle() {
	vec2 pixPos = gl_PointCoord.st * vec2(2.0f, -2.0f) - vec2(1.0f, -1.0f);
	if (length(pixPos) >= 1.0f) discard;
	if (abs(pixPos.x - pixPos.y)>0.2f && abs(pixPos.x + pixPos.y)>0.2f) discard;
}

void paintRGB() {
	float range = sRangeMax - sRangeMin;
	float s_normalized = (fS1 - sRangeMin) / range;
	if(fType == 2) {
		color = vec4(0.1f, 0.1f, 0.1f, 1.0f);
	}
	else if (fType == 3) {
		color = vec4(1.0f, 1.0f, 0.0f, 1.0f);
	}
	else if (fType == 4) {
		color = vec4(0.0f, 1.0f, 1.0f, 1.0f);
	}
	else if (fType == 0) {
		if (s_normalized >= 0.0f && s_normalized < 0.5f)
			color = 2.0f * ((0.5f - s_normalized)* blue + s_normalized* green);
		else if (s_normalized >= 0.5f && s_normalized < 1.0f)
			color = 2.0f * ((1.0f - s_normalized)* green + (s_normalized - 0.5f)* red);
		else if (s_normalized < 0.0f)
			color = blue;
		else
			color = red;
	}
	color.a = 1.0f;
}

void paintRB() {
	float range = sRangeMax - sRangeMin;
	float s_normalized = 2.0f* (fS1 - sRangeMin) / range - 1.0f;
	if (fType == 2) {
		color = vec4(0.1f, 0.1f, 0.1f, 1.0f);
	}
	else if (fType == 3) {
		color = vec4(1.0f, 1.0f, 0.0f, 1.0f);
	}
	else if (fType == 4) {
		color = vec4(0.0f, 1.0f, 1.0f, 1.0f);
	}
	else if (fType == 0) {
		if (s_normalized >= -1.0f && s_normalized < 0.0f)
			color = white - ((0.0f - s_normalized)* (red + green));
		else if (s_normalized >= 0.0f && s_normalized < 1.0f)
			color = white - ((s_normalized - 0.0f)* (green + blue));
		else if (s_normalized < -1.0f)
			color = blue;
		else
			color = red;
	}
	color.a = 1.0f;
}

);

#define GLSL(version, shader)  "#version " #version "\n" #shader

const GLchar* fragment_colorPick = GLSL(330 core,

uniform vec4 colorID;
out vec4 color;

void main() {
	color = colorID;
	return;
}
);