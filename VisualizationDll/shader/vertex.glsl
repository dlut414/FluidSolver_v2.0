/*
 * LICENCE
 * copyright 2014 ~ ****
 * Some rights reserved.
 * Author: HUFANGYUAN
 * Released under CC BY-NC
*/
#version 330 core

layout(location = 0) in int vType;
layout(location = 1) in vec3 vPos;
layout(location = 2) in float v3;
layout(location = 3) in float v4;

flat out int fType;
out vec4 fPos;
out float f3;
out float f4;

uniform mat4 vMvp;

void main()
{
	fType = vType;
    fPos = vMvp * vec4(vPos, 1.0f);
	f3 = v3;
	f4 = v4;

    gl_Position = fPos;
	gl_PointSize = 4.;// +scal / 5.e10;
}
