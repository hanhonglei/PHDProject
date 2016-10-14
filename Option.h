#pragma once

class CSkeOption
{
public:

	bool addNoise;
	bool applyEmbedding;
	bool applyRootFinding;
	bool applySimplification;
	double areaRatioThreshold;
	bool displayIntermediateMesh;
	double laplacianConstraintScale;
	double laplacianConstraintWeight;
	int maxIterations;
	double noiseRatio;
	int numOfImprovement;
	double originalPositionalConstraintWeight;
	double positionalConstraintScale;
	double positionalConstraintWeight;
	bool postSimplify;
	double postSimplifyErrorRatio;
	double shapeEnergyWeight;
	int targetVertexCount;
	bool twoDimensionModel;
	bool useBoundaryVerticesOnly;
	bool useIterativeSolver;
	bool useSamplingEnergy ;
	bool useShapeEnergy;
	bool useSymbolicSolver;
	double volumneRatioThreashold;

	// 保存未收缩的模型体积 [7/5/2011 Han Honglei]
	double originalVolume;

	bool firstTime;
	bool reserved;

	CSkeOption(void);
	~CSkeOption(void);
};
