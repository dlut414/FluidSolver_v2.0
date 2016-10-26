/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
#define GLSL(version, shader)  "#version " #version "\n" #shader

const GLchar* vertex0 = GLSL(330 core,

uniform mat4 vMvp;
uniform mat4 vModelMat;
uniform mat4 vViewMat;
uniform mat4 vProjectionMat;

layout(location = 0) in int vType;
layout(location = 1) in vec3 vPos;
layout(location = 2) in float vS1;
layout(location = 3) in float vS2;

flat out int fType;
out vec4 fPos;
out float fS1;
out float fS2;

void main() {
	vec3 pos = vec3(vPos.x, vPos.y, 0.2f*vS1);
	fType = vType;
	fPos = vMvp * vec4(pos, 1.0f);

	fS1 = vS1;
	fS2 = vS2;

	gl_Position = fPos;
	gl_PointSize = 4.0f;// +scal / 5.e10;
}

);

const GLchar* vertex1 = GLSL(330 core,

uniform mat4 vMvp;
uniform mat4 vModelMat;
uniform mat4 vViewMat;
uniform mat4 vProjectionMat;

layout(location = 0) in int vType;
layout(location = 1) in float vPosX;
layout(location = 2) in float vPosY;
layout(location = 3) in float vS1;
layout(location = 4) in float vS2;

flat out int fType;
out vec4 fPos;
out float fS1;
out float fS2;

void main() {
	vec3 pos = vec3(vPosX, vPosY, 0.0f);
	fType = vType;
	fPos = vMvp * vec4(pos, 1.0f);

	fS1 = vS1;
	fS2 = vS2;

	gl_Position = fPos;
	gl_PointSize = 4.0f;
}

);

#define GLSL(version, shader)  "#version " #version "\n" #shader

const GLchar* vertex_colorPick = GLSL(330 core,

uniform mat4 vMvp;

layout(location = 1) in float vPosX;
layout(location = 2) in float vPosY;

void main() {
	vec3 pos = vec3(vPosX, vPosY, 0.0f);
	vec4 fPos = vMvp * vec4(pos, 1.0f);
	gl_Position = fPos;
	gl_PointSize = 4.0f;
}
);