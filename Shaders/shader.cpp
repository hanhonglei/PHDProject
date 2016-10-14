#include "shader.h"

CGcontext Shader::s_context = NULL;
int Shader::s_shadersLoaded = 0;

Shader::Shader()
:	m_gotFragFile(false)
{
	++s_shadersLoaded;
	m_filenameV[0] = '\0';
	m_filenameF[0] = '\0';
}

Shader::~Shader()
{
	--s_shadersLoaded;
	if(s_shadersLoaded==0 && s_context != NULL)
		cleanUp();
}



void Shader::cleanUp()
{
	cgDestroyContext(s_context);
	s_context = NULL;
	m_filenameV[0] = '\0';
	m_filenameF[0] = '\0';
}


//------------------------------------------------------------------------------------
// Loads and compiles the shader specified. If fragFile is NULL, no fragment shader
// will be used. If both vertex and fragment shaders are in the same file, send that
// filename to both vertFile and fragFile.
void Shader::init(char* vertFile, char* fragFile, char* vertMainFunctionName, char* fragMainFunctionName)
{
	strcpy(m_filenameV, vertFile);
	if(fragFile)
		strcpy(m_filenameF, fragFile);

	//Create context	
	if(s_context == NULL)
		s_context = cgCreateContext();
	if(s_context == NULL)
		MessageBox(NULL, "Can't create shader context", "Error in Shader::init", MB_OK | MB_ICONERROR);


	//Initialize vertex profile
	m_vertProfile = cgGLGetLatestProfile(CG_GL_VERTEX);
	if(m_vertProfile == CG_PROFILE_UNKNOWN)
		MessageBox(NULL, "Can't get vertex shader profile, your graphics card doesn't seem to support vertex shaders in Cg", "Error", MB_OK | MB_ICONERROR);
	cgGLSetOptimalOptions(m_vertProfile);
	
	//Initialize fragment profile
	if(fragFile!=NULL)
	{
		m_fragProfile = cgGLGetLatestProfile(CG_GL_FRAGMENT);
		if(m_fragProfile == CG_PROFILE_UNKNOWN)
			MessageBox(NULL, "Can't get fragment shader profile, your graphics card doesn't seem to support fragment shaders in Cg", "Error", MB_OK | MB_ICONERROR);
		cgGLSetOptimalOptions(m_fragProfile);
	}
	
	
	//Compile and load the vertex program
	compileAndLoad(m_vertProgram, m_vertProfile, vertFile, vertMainFunctionName);
	
	//Compile and load the fragment program
	if(fragFile != NULL)
	{
		compileAndLoad(m_fragProgram, m_fragProfile, fragFile, fragMainFunctionName);
		m_gotFragFile = true;
	}
	
	//test if the profiles are supported on the current system
	if(cgGLIsProfileSupported(m_vertProfile) == CG_FALSE)
		MessageBox(NULL, "Vertex shader profile is not supported on the current system", "WARNING", MB_OK | MB_ICONEXCLAMATION);

	if(cgGLIsProfileSupported(m_fragProfile) == CG_FALSE)
		MessageBox(NULL, "Fragment shader profile is not supported on the current system", "WARNING", MB_OK | MB_ICONEXCLAMATION);
}


//------------------------------------------------------------------------------------
void Shader::compileAndLoad(CGprogram& destinationProgram, CGprofile& profile, char* sourceFile, char* mainFunctionName)
{
	destinationProgram = cgCreateProgramFromFile(s_context, CG_SOURCE, sourceFile, profile, mainFunctionName, 0);
	if(destinationProgram == NULL)	//something went wrong, display error message
	{
		CGerror Error = cgGetError();

		char temptext[100];
		strcpy(temptext, "Error in shader: ");
		strcat(temptext, sourceFile);
		strcat(temptext, "\n\n");
		strcat(temptext, cgGetErrorString(Error));
		MessageBox(NULL, temptext, "Error in shader", MB_OK | MB_ICONERROR);
		
		const char* compileError = cgGetLastListing(s_context);
		MessageBox(NULL, compileError, "Shader compile error", MB_OK);
		return;
	}
	//Load the program
	cgGLLoadProgram(destinationProgram);
}


//------------------------------------------------------------------------------------
// Assign the CGparameter to a handle from the Cg program parameter
void Shader::assignShaderParameter(CGparameter &target, char* name) const
{
	target = cgGetNamedParameter(m_vertProgram, name);
	if (target == NULL)
		target = cgGetNamedParameter(m_fragProgram, name);
	if (target == NULL)
	{
		if(m_filenameV[0] == '\0' && m_filenameF[0] == '\0')
			MessageBox(NULL, "Can't set parameter on uninitialized shader", "Shader::assignShaderParameter, name not found", MB_OK | MB_ICONWARNING);
		else
		{
			char temptext[100];
			strcpy(temptext, "Can't find parameter: ");
			strcat(temptext, name);
			strcat(temptext, "\nIn shader: ");
			strcat(temptext, m_filenameV);
			if(m_filenameF[0] != '\0')
			{
				strcat(temptext, "/");
				strcat(temptext, m_filenameF);
			}
			MessageBox(NULL, temptext, "Shader::assignShaderParameter, name not found", MB_OK | MB_ICONWARNING);
		}
	}
}


//------------------------------------------------------------------------------------
//start using shader
void Shader::begin() const
{
	cgGLEnableProfile(m_vertProfile);
	if(m_gotFragFile)
		cgGLEnableProfile(m_fragProfile);
	
	cgGLBindProgram(m_vertProgram);
	if(m_gotFragFile)
		cgGLBindProgram(m_fragProgram);
}


//------------------------------------------------------------------------------------
//stop using shader
void Shader::end() const
{
	cgGLDisableProfile(m_vertProfile);
	if(m_gotFragFile)
		cgGLDisableProfile(m_fragProfile);
}


//------------------------------------------------------------------------------------
// Use this to send a texture to the Cg-shader, Cg-style (not using glBind).
// Can be set only once because of parameter shadowing
void Shader::setTextureParam(char* name, GLuint textureID) const
{
	CGparameter target;
	assignShaderParameter(target, name);
	cgGLSetTextureParameter(target, textureID);
	cgGLEnableTextureParameter(target);
}


//------------------------------------------------------------------------------------
void Shader::setTextureParam(CGparameter &target, GLuint textureID) const
{
	cgGLSetTextureParameter(target, textureID);
	cgGLEnableTextureParameter(target);
}


//------------------------------------------------------------------------------------
void Shader::disableTextureParam(char* name) const
{
	CGparameter target;
	assignShaderParameter(target, name);
	cgGLDisableTextureParameter(target);
}


//------------------------------------------------------------------------------------
void Shader::disableTextureParam(CGparameter &target) const
{
	cgGLDisableTextureParameter(target);
}


//------------------------------------------------------------------------------------
// Gets the specified matrix from OpenGL and sends it to the shader
// matrix should any of:
//  CG_GL_MODELVIEW_PROJECTION_MATRIX 
//  CG_GL_MODELVIEW_MATRIX
//  CG_GL_PROJECTION_MATRIX
//  CG_GL_TEXTURE_MATRIX 
// transform should any of:
//  CG_GL_MATRIX_IDENTITY
//  CG_GL_MATRIX_TRANSPOSE
//  CG_GL_MATRIX_INVERSE
//  CG_GL_MATRIX_INVERSE_TRANSPOSE 
void Shader::setUniformGLMatrix(char* name, CGGLenum matrix, CGGLenum transform) const
{
	CGparameter target;
	assignShaderParameter(target, name);
	cgGLSetStateMatrixParameter(target,  matrix, transform);
}


//------------------------------------------------------------------------------------
void Shader::setUniformGLMatrix(CGparameter &target, CGGLenum matrix, CGGLenum transform) const
{
	cgGLSetStateMatrixParameter(target,  matrix, transform);
}


//------------------------------------------------------------------------------------
void Shader::setUniform4x4f(char* name, const float* sixteenValues, CGGLenum transform) const
{
	CGparameter target;
	assignShaderParameter(target, name);
	
	glPushMatrix();
		//set matrix
		glMatrixMode(GL_MODELVIEW);
		glLoadMatrixf(sixteenValues);
		
		//send it to the shader
		cgGLSetStateMatrixParameter(target, CG_GL_MODELVIEW_MATRIX, transform);
	glPopMatrix();
}


//------------------------------------------------------------------------------------
void Shader::setUniform4x4f(CGparameter &target, const float* sixteenValues, CGGLenum transform) const
{
	glPushMatrix();
		//set matrix
		glMatrixMode(GL_MODELVIEW);
		glLoadMatrixf(sixteenValues);
		
		//send it to the shader
		cgGLSetStateMatrixParameter(target, CG_GL_MODELVIEW_MATRIX, transform);
	glPopMatrix();
}


//------------------------------------------------------------------------------------
void Shader::setUniform4fv(char* name, const float* fourValues) const
{
	CGparameter target;
	assignShaderParameter(target, name);
	cgGLSetParameter4fv(target, fourValues);
}


//------------------------------------------------------------------------------------
void Shader::setUniform4fv(CGparameter &target, const float* fourValues) const
{
	cgGLSetParameter4fv(target, fourValues);
}


//------------------------------------------------------------------------------------
void Shader::setUniform4f(char* name, const VECTOR4D* fourValues) const
{
	CGparameter target;
	assignShaderParameter(target, name);
	cgGLSetParameter4f(target, fourValues->x, fourValues->y, fourValues->z, fourValues->w);
}


//------------------------------------------------------------------------------------
void Shader::setUniform4f(CGparameter &target, const VECTOR4D* fourValues) const
{
	cgGLSetParameter4f(target, fourValues->x, fourValues->y, fourValues->z, fourValues->w);
}


//------------------------------------------------------------------------------------
void Shader::setUniform3fv(char* name, const float* threeValues) const
{
	CGparameter target;
	assignShaderParameter(target, name);
	cgGLSetParameter3fv(target, threeValues);
}


//------------------------------------------------------------------------------------
void Shader::setUniform3fv(CGparameter &target, const float* threeValues) const
{
	cgGLSetParameter3fv(target, threeValues);
}


//------------------------------------------------------------------------------------
void Shader::setUniform3f(char* name, const VECTOR3D* threeValues) const
{
	CGparameter target;
	assignShaderParameter(target, name);
	cgGLSetParameter3f(target, threeValues->x, threeValues->y, threeValues->z);
}


//------------------------------------------------------------------------------------
void Shader::setUniform3f(CGparameter &target, const VECTOR3D* threeValues) const
{
	cgGLSetParameter3f(target, threeValues->x, threeValues->y, threeValues->z);
}


//------------------------------------------------------------------------------------
void Shader::setUniformInt(char* name, const int* value) const
{
	CGparameter target;
	assignShaderParameter(target, name);
	cgSetParameterValueir(target, 1, value);
}


//------------------------------------------------------------------------------------
void Shader::setUniformInt(CGparameter &target, const int* value) const
{
	cgSetParameterValueir(target, 1, value);
}


//------------------------------------------------------------------------------------
void Shader::setUniformFloat(char* name, const float* value) const
{
	CGparameter target;
	assignShaderParameter(target, name);
	cgSetParameterValuefr(target, 1, value);
}


//------------------------------------------------------------------------------------
void Shader::setUniformFloat(CGparameter &target, const float* value) const
{
	cgSetParameterValuefr(target, 1, value);
}


//------------------------------------------------------------------------------------
void Shader::setUniformToConstant(CGparameter &target)
{
	cgSetParameterVariability(target, CG_LITERAL);
}


//------------------------------------------------------------------------------------
void Shader::setUniformToConstant(char* name)
{
	CGparameter target;
	assignShaderParameter(target, name);
	cgSetParameterVariability(target, CG_LITERAL);
}