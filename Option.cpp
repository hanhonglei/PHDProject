#include "Option.h"

CSkeOption::CSkeOption(void)
{
	addNoise;
	applyEmbedding = true;
	applyRootFinding;
	applySimplification = true;
	areaRatioThreshold = 0.001;
	displayIntermediateMesh;
	laplacianConstraintScale = 2.0;
	laplacianConstraintWeight = 1.0;
	maxIterations = 30;
	noiseRatio = 0.02;
	numOfImprovement = 100;
	originalPositionalConstraintWeight;
	positionalConstraintScale = 1.5;
	positionalConstraintWeight = 1.0;
	postSimplify = true;
	postSimplifyErrorRatio = 0.9;
	shapeEnergyWeight = 0.1;
	targetVertexCount = 10;
	twoDimensionModel;
	useBoundaryVerticesOnly = true;
	useIterativeSolver;
	useSamplingEnergy = true;
	useShapeEnergy = true;
	useSymbolicSolver;
	volumneRatioThreashold = 1E-05;
	originalVolume = 0;
	firstTime = true;
}

CSkeOption::~CSkeOption(void)
{
}
