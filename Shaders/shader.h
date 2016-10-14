/**************************************************************************
Shader
	Class to use when using Cg-shaders. Init with [init], then use
	[begin]/[end] when using the shader on an object.
		If any parameters are to be sent to the shader then use
	the setUniform-functions which take either a CGparameter or
	a char*. By only using the functions taking char*, Cg is
	completely encapsulated. If, on the other hand, you want to
	use the CGparameters directly, this is possible too, but then
	you need to assign them first using [assignShaderParameter] and
	[assignShaderTextureParameter]. This method might be slightly
	faster since [assignShader_] only has to be called once.
	

	Copyright 2007, Anders Nivfors.
	
	Permission is granted to anyone to use this software for any purpose,
	including commercial applications, and to alter it and redistribute it
	freely, as long as the above copyright notice and this disclaimer are
	retained. Furthermore, if you use this software in a product, an
	acknowledgment in the product documentation would be appreciated but is
	not required.
***************************************************************************/

#ifndef SHADER_H
#define SHADER_H


#include <windows.h>
#include "maths.h"
#include <cg/cg.h>
#include <cg/cgGL.h>


class Shader
{
	public:
		Shader();
		~Shader();
		void init(char* vertFile, char* fragFile=NULL, char* vertMainFunctionName="mainV", char* fragMainFunctionName="mainF");
		void assignShaderParameter(CGparameter &target, char* name) const;			//assign a parameter (or use the setUniform-functions instead)
		void begin() const;	//start using shader
		void end() const;	//stop using shader
		void setTextureParam(char* name, GLuint textureID) const;					//"bind" a texture (instead of using TEXUNITx)
		void setTextureParam(CGparameter &target, GLuint textureID) const;			//"bind" a texture (instead of using TEXUNITx)
		void disableTextureParam(char* name) const;									//"unbind" a texture
		void disableTextureParam(CGparameter &target) const;						//"unbind" a texture
		void setUniformGLMatrix(char* name, CGGLenum matrix, CGGLenum transform=CG_GL_MATRIX_IDENTITY) const;			//send a specified openGL matrix to the shader
		void setUniformGLMatrix(CGparameter &target, CGGLenum matrix, CGGLenum transform=CG_GL_MATRIX_IDENTITY) const;	//send a specified openGL matrix to the shader
		void setUniform4x4f(char* name, const float* sixteenValues, CGGLenum transform=CG_GL_MATRIX_IDENTITY) const;	//send a custom 4x4 matrix to the shader
		void setUniform4x4f(CGparameter &target, const float* sixteenValues, CGGLenum transform=CG_GL_MATRIX_IDENTITY) const;	//send a custom 4x4 matrix to the shader
		void setUniform4fv(char* name, const float* fourValues) const;				//send 4 floats to the shader
		void setUniform4fv(CGparameter &target, const float* fourValues) const;		//send 4 floats to the shader
		void setUniform4f(char* name, const VECTOR4D* fourValues) const;			//send 4 floats to the shader (using a VECTOR4D)
		void setUniform4f(CGparameter &target, const VECTOR4D* fourValues) const;	//send 4 floats to the shader (using a VECTOR4D)
		void setUniform3fv(char* name, const float* threeValues) const;				//send 3 floats to the shader
		void setUniform3fv(CGparameter &target, const float* threeValues) const;	//send 3 floats to the shader
		void setUniform3f(char* name, const VECTOR3D* threeValues) const;			//send 3 floats to the shader (using a VECTOR3D)
		void setUniform3f(CGparameter &target, const VECTOR3D* threeValues) const;	//send 3 floats to the shader (using a VECTOR3D)
		void setUniformInt(char* name, const int* value) const;						//send an int to the shader
		void setUniformInt(CGparameter &target, const int* value) const;			//send an int to the shader
		void setUniformFloat(char* name, const float* value) const;					//send a float to the shader
		void setUniformFloat(CGparameter &target, const float* value) const;		//send a float to the shader
		void setUniformToConstant(CGparameter &target);								//set the param to be a constant (program needs to be compiled if param changes)
		void setUniformToConstant(char* name);										//set the param to be a constant (program needs to be compiled if param changes)
		
	protected:
		void cleanUp();
		void compileAndLoad(CGprogram& destinationProgram, CGprofile& profile, char* sourceFile, char* mainFunctionName);
		CGprogram			m_vertProgram;
		CGprogram			m_fragProgram;
		CGprofile			m_vertProfile;
		CGprofile			m_fragProfile;
		bool				m_gotFragFile;		//true if a fragment shader was loaded in init
		char				m_filenameV[50];	//only for debug purpose (used in error messages)
		char				m_filenameF[50];	//only for debug purpose (used in error messages)
		//static variables
		static CGcontext	s_context;			//Cg context (container for multiple Cg programs)
		static int			s_shadersLoaded;	//keeps count on how many shaders we have created, when the last is deleted, the context is deleted
};


#endif