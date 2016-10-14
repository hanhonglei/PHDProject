// SegmentParam.cpp : implementation file
//

#include "stdafx.h"
#include "MeshSimp.h"
#include "SegmentParam.h"
#include "afxdialogex.h"

// SegmentParam dialog

IMPLEMENT_DYNAMIC(SegmentParam, CDialogEx)

SegmentParam::SegmentParam(CWnd* pParent /*=NULL*/)
	: CDialogEx(SegmentParam::IDD, pParent)
{

}

SegmentParam::~SegmentParam()
{
}

void SegmentParam::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(SegmentParam, CDialogEx)
	ON_EN_CHANGE(IDC_SIGMA, &SegmentParam::OnEnChangeSigma)
	ON_WM_DESTROY()
END_MESSAGE_MAP()


// SegmentParam message handlers


void SegmentParam::OnEnChangeSigma()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialogEx::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
}


BOOL SegmentParam::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// TODO:  Add extra initialization here
	CString   str; 
	str.Format("%f", m_sigma); 
	SetDlgItemText(IDC_SIGMA, str);

	str.Format("%f", m_k); 
	SetDlgItemText(IDC_K, str);

	str.Format("%f", m_min_size); 
	SetDlgItemText(IDC_MIN, str);

	str.Format("%d", m_min_size); 
	CheckDlgButton(IDC_CHECK1, m_bFace);

	return TRUE;  // return TRUE unless you set the focus to a control
	// EXCEPTION: OCX Property Pages should return FALSE
}

void SegmentParam::SetValues(float sigma, float k, float min_size, bool bFace)
{
	m_sigma = sigma;
	m_k = k;
	m_min_size = min_size;
	m_bFace = bFace;
}
void SegmentParam::GetValues(float &sigma, float &k, float &min_size, bool &bFace)
{
	sigma = m_sigma;
	k = m_k;
	min_size = m_min_size;	
	bFace = m_bFace;
}


void SegmentParam::OnDestroy()
{
	CDialogEx::OnDestroy();

	// TODO: Add your message handler code here
	CString   str; 

	GetDlgItemText(IDC_SIGMA, str);
	m_sigma =  (float)atof(str.GetBuffer(str.GetLength()));

	GetDlgItemText(IDC_K, str);
	m_k =  (float)atof(str.GetBuffer(str.GetLength()));

	GetDlgItemText(IDC_MIN, str);
	m_min_size =  (float)atof(str.GetBuffer(str.GetLength()));

	m_bFace = IsDlgButtonChecked(IDC_CHECK1);

}
