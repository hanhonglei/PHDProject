// MeshSimp.h : main header file for the MeshSimp application
//
#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"       // main symbols


// CMeshSimpApp:
// See MeshSimp.cpp for the implementation of this class
//

class CMeshSimpApp : public CWinApp
{
public:
	CMeshSimpApp();


// Overrides
public:
	virtual BOOL InitInstance();

// Implementation
	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()
	virtual BOOL OnCmdMsg(UINT nID, int nCode, void* pExtra, AFX_CMDHANDLERINFO* pHandlerInfo);
};

extern CMeshSimpApp theApp;