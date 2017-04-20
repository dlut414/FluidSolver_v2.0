/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
#pragma once
#include "Header.h"
#include "Camera.h"

namespace VIS {

	enum DISPLAYMODE { DMODE_ONE = 1, DMODE_TWO = 2, DMODE_THREE = 3, DMODE_FOUR = 4, DMODE_FIVE = 5, DMODE_SIX = 6, };

	class Controller {
	public:
		Controller() {
			m_mode = DMODE_ONE;
			f_scaleVel = 0.0010f;
			f_panVel = 0.01f;
			f_visScaleVel = 0.1f;

			m_scale = m_initScale = glm::vec3(1.f);
			m_rotation = m_initRotation = glm::angleAxis<float>(-glm::pi<float>() * 0.f, glm::vec3(1, 0, 0));;
			m_pan = m_initPan = glm::vec3(0.f, 0.f, 0.f);

			m_initCameraPosition = glm::vec3(0.0f, 0.0f, 1.0f);
			m_initCameraRotation = glm::angleAxis<float>(glm::pi<float>() * 0.0f, glm::vec3(1, 0, 0));

			f_pointRadius = 2.0f;
			f_pointScale = 1.0f;
			f_near = 0.001f;
			f_far = 1000.f;

			i_init = 0;
			i_dirty = 1;
			i_stop = 1;
			i_leave = 0;
			i_point = 0;
			i_save = 0;
			i_sens = 0;
			i_bmp = 0;
			i_senSwitch = 1;
			i_bmpSwitch = 1;
			u_width = 1800;
			u_height = 250;
			u_width_max = 1800;
			u_height_max = 1000;
			f_sRangeMax = 1.0f;
			f_sRangeMin = -1.0f;

			m_camera.SetPosition(m_initCameraPosition);
			//m_camera.SetProjectionRH(45.0f, float(u_width)/float(u_height), f_near, f_far);
			m_camera.SetProjectionOR(-0.55f*float(u_width) / float(u_height), 0.55f*float(u_width) / float(u_height), -0.55f, 0.55f, f_near, f_far);
		}
		~Controller() {}

		void setProjectionOR(float left, float right, float bottom, float top) {
			const float whRatio = float(u_width) / float(u_height);
			const float centerX = (left + right) / 2;
			const float dX = whRatio* (top - bottom);
			m_camera.SetProjectionOR(centerX - dX / 2, centerX + dX / 2, bottom, top, f_near, f_far);
		}

		void clickMouse(int button, int state, int x, int y) {
			m_mousePos = glm::ivec2(x, y);
			if (state == GLUT_UP) return;
			switch (button) {
			case GLUT_LEFT_BUTTON:
			{
				i_mouseButton = GLUT_LEFT_BUTTON;
				break;
			}
			case GLUT_RIGHT_BUTTON:
			{
				i_mouseButton = GLUT_RIGHT_BUTTON;
				break;
			}
			case GLUT_MIDDLE_BUTTON:
			{
				i_mouseButton = GLUT_MIDDLE_BUTTON;
				break;
			}
			}
		}
		void moveMouse(int x, int y) {
			glm::ivec2 mousePos = glm::ivec2(x, y);
			glm::vec2 delta = glm::vec2(mousePos - m_mousePos);
			m_mousePos = mousePos;

			switch (i_mouseButton) {
			case GLUT_LEFT_BUTTON:
			{
				glm::quat rotX = glm::angleAxis<float>(glm::radians(delta.y) * 0.5f, glm::vec3(1, 0, 0));
				glm::quat rotY = glm::angleAxis<float>(glm::radians(delta.x) * 0.5f, glm::vec3(0, 1, 0));
				m_rotation = (rotX * rotY) * m_rotation;
				break;
			}
			case GLUT_RIGHT_BUTTON:
			{
				m_pan += glm::vec3(f_panVel*delta.x, -f_panVel*delta.y, 0.0f);
				break;
			}
			case GLUT_MIDDLE_BUTTON:
			{
				m_scale += glm::vec3(delta.y * f_scaleVel);
				m_scale = glm::max(m_scale, glm::vec3(0.f, 0.f, 0.f));
				break;
			}
			}
			i_dirty = 1;
		}
		void rollMouse(int button, int dir, int x, int y) {
			m_scale *= dir * f_scaleVel;
			i_dirty = 1;
		}

		void reshapeWindow() {
			i_dirty = 1;
		}

		void pressKey(unsigned char key, int a, int b) {
			switch (key) {
			case 0x1b: //esc
				i_leave = 1;
				break;
			case 0x0d: //enter
				i_stop = !i_stop;
				break;
			case 0x70: //p
				i_point = !i_point;
				break;
			case 0x20: //space
				m_scale = m_initScale;
				m_rotation = m_initRotation;
				m_pan = m_initPan;
				break;
			case 0x2c: //,
				i_init = 1;
				break;
			case 0x2e: //.
				i_init = 1;
				break;
			case 0x53: //S
				i_save = 1;
				break;
			case 0x73: //s
				i_sens = !i_sens;
				break;
			case 0x10: //ctrl p
				i_bmp = 1;
				break;
			default:
				break;
			}
			i_dirty = 1;
		}

	public:
		DISPLAYMODE m_mode;
		Camera      m_camera;
		glm::ivec2  m_mousePos;
		glm::quat   m_rotation;
		glm::vec3   m_scale;
		glm::vec3   m_pan;
		glm::mat4   m_mvp;
		glm::mat4   m_mvpInv;
		glm::mat4   m_modelMat;
		glm::mat4   m_viewMat;
		glm::mat4   m_projectionMat;
		glm::mat4   m_viewModelMat;
		glm::mat4   m_projectionMatInv;

		GLfloat     f_pointRadius;
		GLfloat     f_pointScale;
		GLfloat     f_near;
		GLfloat     f_far;
		GLuint      u_width;
		GLuint      u_height;
		GLuint		u_width_max;
		GLuint		u_height_max;

		int			i_init;
		int			i_dirty;
		int			i_stop;
		int			i_leave;
		int			i_point;
		int			i_save;
		int			i_sens;
		int			i_bmp;
		int			i_senSwitch;
		int			i_bmpSwitch;
		int         i_mouseButton;
		int         i_file;

		float		f_sRangeMax;
		float		f_sRangeMin;

	private:
		float       f_scaleVel;
		float       f_panVel;
		glm::vec3   m_initCameraPosition;
		glm::quat   m_initCameraRotation;
		glm::quat   m_initRotation;
		glm::vec3   m_initScale;
		glm::vec3   m_initPan;
		float		f_visScaleVel;
	};

}

