#ifndef HEADER_PCH_H_INCLUDED
#define HEADER_PCH_H_INCLUDED

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

//#define GLEW_STATIC
#include <GL/glew.h>

#ifdef _WIN32
#include <GL/wglew.h> // For wglSwapInterval
#endif

//#define FREEGLUT_STATIC
#include <GL/freeglut.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/quaternion.hpp>

#include <AntTweakBar.h>

#define __FILENAME__ (strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\') + 1 : __FILE__)
#define PRINT(x) std::cout << __FILENAME__ << " @Ln " << __LINE__ << ": " << #x << " = " << x << std::endl;

#endif // HEADER_PCH_H_INCLUDED
