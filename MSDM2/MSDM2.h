#ifndef HEADER_MSDM2_H
#define HEADER_MSDM2_H
#include "..\jmspmesh\vec3.h"
#include "stdmix.h"
#include <MxQSlim.h>
#include <MxSMF.h>

class MSDM2 
{
	//Q_OBJECT
	//Q_INTERFACES(mepp_component_plugin_interface);

	public:
		MSDM2();
		~MSDM2() {};

	public :
		static double MSDM2_computation(MxStdModel *pSrcMesh, MxStdModel *pSimpMesh, bool bSym = false);
		static void DistanceToColorMap();


		
		/**
		 \fn	
		
		 \brief	Computes the multiscale MSDM2 metric
				
		 \param	pSrcMesh	  	The source Polyhedron.
		 \param	pDesMesh	  	The Polyhedron to be calced.
		 \param	NbLevel		Number of scales used
		 \param	maxdim	The max dimension of the Bounding Box Length.
		 \param [out] MSDM2Value	The computed value of the MSDM2 metric
		
		
		 */
		static double ProcessMSDM2_Multires(MxStdModel *pSrcMesh, MxStdModel *pDesMesh,int NbLevel, double maxdim,double & FastMSDM);
		
		/**
		 \fn	void MSDM2_Component::Matching_Multires_Init(PolyhedronPtr m_PolyDegrad , PolyhedronPtr m_PolyOriginal, 
				double Length,Facet * _TabMatchedFacet);
				
		
		 \brief	Initialize the matching process
				
		 \param	m_PolyDegrad	  	The first Polyhedron.
		 \param	m_PolyOriginal	  	The second Polyhedron.
		 \param	[out] _TabMatchedFacet	Facets from m_PolyOriginal on which vertices from m_PolyDegrad are projected
				
		
		 */		 
		static void Matching_Multires_Init( MxStdModel *pDesMesh, MxStdModel *pSrcMesh/*, MxFace * _TabMatchedFacet*/);
		
		/**
		 \fn	void MSDM2_Component::Matching_Multires_Update(PolyhedronPtr m_PolyDegrad , PolyhedronPtr m_PolyOriginal, 
				double Length,Facet * _TabMatchedFacet);
				
		
		 \brief Updates the matching process
				
		 \param	m_PolyDegrad	  	The first Polyhedron.
		 \param	_TabMatchedFacet	Facets from m_PolyOriginal on which vertices from m_PolyDegrad are projected
				
		
		 */	
		static void Matching_Multires_Update( MxStdModel *pDesMesh, MxStdModel *pSrcMesh, int nScale = -1);
		static void Matching_Multires_Update( MxStdModel *pDesMesh, MxStdModel *pSrcMesh, int nScale, int vID);

		/**
		 \fn	void MSDM_Component::ProcessMSDM2_per_vertex(Vertex_iterator pVertex,double radius,
				std::vector<double> & TabDistance1,std::vector<double>& TabDistance2,std::vector<Point3d> &TabPoint1,std::vector<Point3d> &TabPoint2);
	
		 \brief	Computes the local neighborhoods
		
		
		 \param	pVertex	The considered vertex
		 \param radius : radius of the neighborhood
		 \param [out] TabDistance1 : Curvature values from the neighborhood of pVertex regarding the first polyhedron
		 \param [out] TabDistance2 : Curvature values from the neighborhood of pVertex regarding the second polyhedron
		 \param [out] TabPoint1 : 3D points from the neighborhoodof pVertex regarding the first polyhedron
		 \param [out] TabPoint2 : 3D points from the neighborhood of pVertex regarding the second polyhedron
		 
		 */
		static void ProcessMSDM2_per_vertex(MxStdModel *pMesh, MxVertexID verID, double radius,std::vector<double> & TabDistance1
			,std::vector<double>& TabDistance2,std::vector<vec3> &TabPoint1,std::vector<vec3> &TabPoint2, int nScale = -1);
	
		/**
		 \fn	double MSDM_Component:: ComputeStatistics(Vertex* pVertex, double Param,std::vector<double> & TabDistance1,
				 std::vector<double>& TabDistance2,std::vector<Point3d> &TabPoint1,std::vector<Point3d> &TabPoint2,
				 double radius,double dim);
		
		 \brief	Calculates the local curvature statistics per vertex.
		
		
		 \param	pVertex	The considered vertex
		
		 \param TabDistance1 : Curvature values from the neighborhood of pVertex regarding the first polyhedron
		 \param TabDistance2 : Curvature values from the neighborhood of pVertex regarding the second polyhedron
		 \param TabPoint1 : 3D points from the neighborhoodof pVertex regarding the first polyhedron
		 \param TabPoint2 : 3D points from the neighborhood of pVertex regarding the second polyhedron
		 \param	radius				   	The radius of the neighborhood.
		 \param	dim	The max dimension of the Bounding Box Length.
		
		
		 */
		static void ComputeStatistics(MxStdModel *pMesh, MxVertexID verID,  double Param,
			std::vector<double> & TabDistance1,std::vector<double>& TabDistance2,std::vector<vec3> &TabPoint1,std::vector<vec3> &TabPoint2,double radius, double dim, int nScale = -1);
		
		/**
		 \fn	void MSDM_Component::ComputeMaxMin(PolyhedronPtr polyhedron_ptr);
		
		 \brief	Calculates the maximum and minimum local MSDM values for rendering MSDM color map
		
		
		 \param	polyhedron_ptr	The polyhedron.
		 \param MetricOrHausdorff : 0 Hausdorff scalar field , 1 MSDM2 scalar field
		 */
		//void ComputeMaxMin(PolyhedronPtr P, int MetricOrHausdorff);
		
		/*!
		* \brief This method map the local MSDM2 scalar field into vertex colors
		* \param pMesh : input polyhedra
		* \param MetricOrHausdorff : 0 Hausdorff scalar field , 1 MSDM2 scalar field
		*/	
		//void ConstructColorMap(PolyhedronPtr P, int MetricOrHausdorff);
		
		/**
		 \fn	void MSDM2_Component::KmaxKmean(PolyhedronPtr polyhedron_ptr,double coef);
		
		 \brief	Computes the mean curvature field (kmin+kmax)/2 and normalize it according to the size of the model
				
		 \param	polyhedron_ptr	The polyhedron.
		 \param	coef		  	The normalization coef.
		 */
		static void KmaxKmean(MxStdModel *pMesh, double coef, int nScale = -1);
		static void KmaxKmean(MxStdModel *pMesh, double coef, int nScale, MxVertexID nV);
		
		/**
		 \fn	double MSDM2_Component::getMaxDim(PolyhedronPtr polyhedron_ptr);
		
		 \brief	Gets the maximum dimension of the object bounding box length.
				
		 \param	polyhedron_ptr	The polyhedron.
		
		 \return	The maximum dimension.
		 */
		static double getMaxDim(MxStdModel *pMesh);

		
		/*! \brief Sets of 3D points of local windows*/
		static std::vector<vec3> TabPoint1;
		static std::vector<vec3> TabPoint2;

		// han
		static void MSDM2SimpInit(MxStdModel *pSrcMesh, MxStdModel *pSimpMesh);
		static void MatchingSimpInit(MxStdModel *pMesh);
		static void UpdateMSDM2(MxStdModel *pSrcMesh, MxStdModel *pSimpMesh, MxVertexID fromVID, MxVertexID toVID, double & FastMSDM ,bool bGlobalMsdm = false, bool bSym = false);
		static void Matching_Multires_Init( MxStdModel *pDesMesh, MxStdModel *pSrcMesh, MxVertexID fromVID, MxVertexID toVID);
		static void GetMSDMInfo(MxStdModel* pMesh, double & FastMSDM ,double &MinLocalMSDM, double &MaxLocalMSDM);
		static double mini_radius;
		static double radius_step;

		private:
		// 更新最近点
		static void FindNearest( MxStdModel *pDesMesh, MxStdModel *pSrcMesh, MxVertexID vID, MxVertexID nearestVID);


};

#endif