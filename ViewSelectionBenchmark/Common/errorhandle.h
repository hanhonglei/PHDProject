#pragma once

/************************************************************************/
/*      Define Cg, OpenGL, FBO error checking handler                   */
/************************************************************************/
bool check_frame_buffer_object_status(const char *label = 0);
bool check_OpenGL_errors(const char *label);