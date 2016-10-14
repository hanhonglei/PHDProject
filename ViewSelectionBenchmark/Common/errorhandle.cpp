#include <gl/glew.h>
#include <cg/cg.h>
#include <cg/cgGL.h>

#include <iostream>
#include "errorhandle.h"


bool check_frame_buffer_object_status(const char* label)
{
    GLenum status;
    status=(GLenum)glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    
    if(label) std::cout<<label<<": ";
    
    switch(status) 
    {
    case GL_FRAMEBUFFER_COMPLETE_EXT:
        std::cout<<"Framebuffer complete\n";
        return true;
    case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:
        std::cout<<"Framebuffer incomplete,incomplete attachment\n";
        return false;
    case GL_FRAMEBUFFER_UNSUPPORTED_EXT:
        std::cout<<"Unsupported framebuffer format\n";
        return false;
    case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
        std::cout<<"Framebuffer incomplete,missing attachment\n";
        return false;
    case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
        std::cout<<"Framebuffer incomplete,attached images must have same dimensions\n";
        return false;
    case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
        std::cout<<"Framebuffer incomplete,attached images must have same format\n";
        return false;
    case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
        std::cout<<"Framebuffer incomplete,missing draw buffer\n";
        return false;
    case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
        std::cout<<"Framebuffer incomplete,missing read buffer\n";
        return false;
    }
    return true;
}

bool check_OpenGL_errors(const char *label) 
{
    GLenum errCode;
    const GLubyte *errStr;
    if ((errCode = glGetError()) != GL_NO_ERROR) {
        errStr = gluErrorString(errCode);
        std::cout<<"OpenGL ERROR:\t"<<errStr<<"\nLabel:\t"<<label<<"\n";
        return false;
    }
    return true;
}