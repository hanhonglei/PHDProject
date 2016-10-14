#pragma once


// SegmentParam dialog

class SegmentParam : public CDialogEx
{
	DECLARE_DYNAMIC(SegmentParam)

public:
	SegmentParam(CWnd* pParent = NULL);   // standard constructor
	virtual ~SegmentParam();

// Dialog Data
	enum { IDD = IDD_DIALOG1 };

protected:
	float m_sigma;
	float m_k;
	float m_min_size;
	bool m_bFace;
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnEnChangeSigma();
	virtual BOOL OnInitDialog();
	void SetValues(float sigma, float k, float min_size, bool bFace);
	void GetValues(float &sigma, float &k, float &min_size, bool &bFace);
	afx_msg void OnDestroy();
};
