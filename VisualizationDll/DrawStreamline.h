/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
#pragma once
#include "Header.h"
#include "Draw_.h"
#include "Shader.h"
#include <Controller.h>
#include <BBox.h>
#define EPS 0.00001

namespace VIS {

	template <typename R>
	class DrawStreamline : public Draw_ {
		template <typename R_ = R>	struct DataType {};
		template <>					struct DataType<float> { enum { value = GL_FLOAT, }; };
		template <>					struct DataType<double> { enum { value = GL_DOUBLE, }; };
	public:
		DrawStreamline() : Draw_() { init(); }
		~DrawStreamline() {}

		template <typename V>
		void makeStreamline(const Controller* const controlPtr, void(*interp)(const V, const V, V&, V&)) {
			//V px, py, qx, qy;
			//int nLine;
			//V res, sLength;
			//std::cout << " Specify 2 points of endpoints of a integration line (px, py, qx, qy): " << std::endl;
			//std::cin >> px >> py >> qx >> qy;
			//std::cout << " Specify the number of streamlines (nLine): " << std::endl;
			//std::cin >> nLine;
			//std::cout << " Specify resolution/length of lines (res): " << std::endl;
			//std::cin >> res;
			//std::cout << " Specify length of streamlines: " << std::endl;
			//std::cin >> sLength;
			linex.clear();
			liney.clear();
			V dx = (controlPtr->v_p2.x - controlPtr->v_p1.x) / (controlPtr->i_nLines + 1);
			V dy = (controlPtr->v_p2.y - controlPtr->v_p1.y) / (controlPtr->i_nLines + 1);
			for (int i = 1; i <= controlPtr->i_nLines; i++) {
				V startx = controlPtr->v_p1.x + dx*i, starty = controlPtr->v_p1.y + dy*i;
				V forwardVel[2];
				interp(startx, starty, forwardVel[0], forwardVel[1]);
				V backwardVel[2] = { forwardVel[0], forwardVel[1] };
				V velmag = sqrt(forwardVel[0] * forwardVel[0] + forwardVel[1] * forwardVel[1]);
				V length = 0;
				V forwardx = startx, forwardy = starty;
				V backwardx = startx, backwardy = starty;
				std::vector<V> forwardLines_x; forwardLines_x.push_back(startx);
				std::vector<V> forwardLines_y; forwardLines_y.push_back(starty);
				linex.push_back(forwardLines_x);
				liney.push_back(forwardLines_y);
				while (velmag > EPS && length <= controlPtr->f_sLength) {
					forwardx += controlPtr->f_dStep*forwardVel[0];
					forwardy += controlPtr->f_dStep*forwardVel[1];
					linex[2*i - 2].push_back(forwardx);
					liney[2*i - 2].push_back(forwardy);
					velmag = sqrt(forwardVel[0] * forwardVel[0] + forwardVel[1] * forwardVel[1]);
					length += controlPtr->f_dStep * velmag;
					interp(forwardx, forwardy, forwardVel[0], forwardVel[1]);
				}
				std::vector<V> backwardLines_x; backwardLines_x.push_back(startx);
				std::vector<V> backwardLines_y; backwardLines_y.push_back(starty);
				linex.push_back(backwardLines_x);
				liney.push_back(backwardLines_y);
				length = 0;
				velmag = sqrt(backwardVel[0] * backwardVel[0] + backwardVel[1] * backwardVel[1]);
				while (velmag > EPS && length <= controlPtr->f_sLength) {
					backwardx -= controlPtr->f_dStep*backwardVel[0];
					backwardy -= controlPtr->f_dStep*backwardVel[1];
					linex[2*i - 1].push_back(backwardx);
					liney[2*i - 1].push_back(backwardy);
					velmag = sqrt(backwardVel[0] * backwardVel[0] + backwardVel[1] * backwardVel[1]);
					length += controlPtr->f_dStep * velmag;
					interp(backwardx, backwardy, backwardVel[0], backwardVel[1]);
				}
			}
		}
		void Draw(const Controller* const controlPtr) const {
			///clear framebuffer
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
			if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) printf("fbo_def not ready\n");
			glViewport(0, 0, controlPtr->u_width, controlPtr->u_height);
			glFrontFace(GL_CCW);
			glClearDepth(1);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			///use program1
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
			if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) printf("0 not ready\n");
			glViewport(0, 0, controlPtr->u_width, controlPtr->u_height);
			glFrontFace(GL_CCW);
			glDepthFunc(GL_LESS);

			glClearDepth(1);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glUseProgram(shaderObj.programID[0]);

			glUniform1f(glGetUniformLocation(shaderObj.programID[0], "sRangeMax"), controlPtr->f_sRangeMax);
			glUniform1f(glGetUniformLocation(shaderObj.programID[0], "sRangeMin"), controlPtr->f_sRangeMin);
			glUniformMatrix4fv(glGetUniformLocation(shaderObj.programID[0], "vMvp"), 1, GL_FALSE, &(controlPtr->m_mvp[0][0]));
			glUniformMatrix4fv(glGetUniformLocation(shaderObj.programID[0], "fMvpInv"), 1, GL_FALSE, &(controlPtr->m_mvpInv[0][0]));
			glUniformMatrix4fv(glGetUniformLocation(shaderObj.programID[0], "vModelMat"), 1, GL_FALSE, &(controlPtr->m_modelMat[0][0]));
			glUniformMatrix4fv(glGetUniformLocation(shaderObj.programID[0], "vViewMat"), 1, GL_FALSE, &(controlPtr->m_viewMat[0][0]));
			glUniformMatrix4fv(glGetUniformLocation(shaderObj.programID[0], "vProjectionMat"), 1, GL_FALSE, &(controlPtr->m_projectionMat[0][0]));

			for (int i = 0; i < linex.size() && i < liney.size(); i++) {
				glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
				glBufferData(GL_ARRAY_BUFFER, linex[i].size()*sizeof(R), linex[i].data(), GL_STATIC_DRAW);
				glVertexAttribPointer(0, 1, DataType<>::value, GL_FALSE, 0, (void*)0);

				glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
				glBufferData(GL_ARRAY_BUFFER, linex[i].size()*sizeof(R), liney[i].data(), GL_STATIC_DRAW);
				glVertexAttribPointer(1, 1, DataType<>::value, GL_FALSE, 0, (void*)0);

				glEnableVertexAttribArray(0);
				glEnableVertexAttribArray(1);

				glDrawArrays(GL_LINE_STRIP, 0, (GLsizei)(linex[i].size()));

				glDisableVertexAttribArray(0);
				glDisableVertexAttribArray(1);
			}

			glUseProgram(0);
		}

	public:
		GLuint bgTex;
		GLuint frontFaceTex;
		GLuint dirTex;
		GLuint dbo[2];

		GLuint volNormTex;
		GLuint volDataTex;
		GLuint cellInFluidTex;

		Shader shaderObj;

		BBox<float> box;
		GLuint ix, iy, iz;
		GLuint icx, icy, icz;

		void* volNorm;
		void* volData;
		void* cellInFluid;

		std::vector<std::vector<R>> linex;
		std::vector<std::vector<R>> liney;

	private:
		void init() {
			Draw_::initVAO(Draw_::vao);
			Draw_::initVBO(Draw_::vbo[0]);
			Draw_::initVBO(Draw_::vbo[1]);
			Draw_::initVBO(Draw_::vbo[2]);
			Draw_::initVBO(Draw_::vbo[3]);
			Draw_::initVBO(Draw_::vbo[4]);
			Draw_::initVBO(Draw_::vbo[5]);
			initShader();
			linex.clear();
			liney.clear();
		}
		void initShader() {
			shaderObj.programID.push_back(shaderObj.LoadShader_inline(vertex_stream, fragment_stream));
		}
	};

}
